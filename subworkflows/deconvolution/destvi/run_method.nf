process buildDestVIModel {
    tag "destvi_build_$tag_suffix"
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    container 'csangara/sp_destvi:latest'
    echo true

    input:
        path (sc_input)
        tuple path (sp_input_h5ad), path (sp_input_rds)
    output:
        tuple path ("adata.h5ad"), path ("model.pt")

    script:
        tag_suffix = file(sc_input).getSimpleName()
        epochs = ( params.epoch_build ==~ /default/ ? "" : "-e $params.epoch_build")
        args = ( params.deconv_args.destvi.build ? params.deconv_args.destvi.build : "" )
        cuda_device = ( params.gpu ? params.cuda_device : "cpu" )
        println ("Building DestVI model with ${ (params.gpu) ? "GPU" : "CPU" }...")
        """
        source activate scvi-env
        python $params.rootdir/subworkflows/deconvolution/destvi/build_model.py \
            $sc_input $sp_input_h5ad $cuda_device -a $params.annot \
            -o \$PWD/model $epochs $args
        mv model/* .
        """

}

process fitDestVIModel {
    tag "destvi_$sp_file_basename"
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    container 'csangara/sp_destvi:latest'
    echo true

    input:
        tuple path (sp_input), path (sp_input_rds)
        tuple path (h5ad_file), path (model_file)
    output:
        tuple val('destvi'), path("$output"), path (sp_input_rds)
    script:
        sp_file_basename = file(sp_input).getSimpleName()
        output = "proportions_destvi_${sp_file_basename}${params.runID_props}.preformat"
        epochs = ( params.epoch_fit ==~ /default/ ? "" : "-e $params.epoch_fit")
        args = ( params.deconv_args.destvi.fit ? params.deconv_args.destvi.fit : "" )
        cuda_device = ( params.gpu ? params.cuda_device : "cpu" )

        println ("Received files $h5ad_file and $model_file")
        println ("Fitting DestVI model with ${ (params.gpu) ? "GPU" : "CPU" }...")
        
        """
        source activate scvi-env
        mkdir model; mv -t model/ $h5ad_file $model_file
        python $params.rootdir/subworkflows/deconvolution/destvi/fit_model.py \
            $sp_input $cuda_device $epochs $args -o \$PWD 
        mv proportions.tsv $output
        """
}
