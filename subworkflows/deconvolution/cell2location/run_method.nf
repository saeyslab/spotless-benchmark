process buildCell2locationModel {
    tag "c2l_build_$tag_suffix"
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    container 'csangara/sp_cell2location:latest'
    echo true

    input:
        path (sc_input)
    output:
        path "sc.h5ad"

    script:
        tag_suffix = file(sc_input).getSimpleName()
        sample_id_arg = ( params.sampleID ==~ /none/ ? "" : "-s $params.sampleID" )
        epochs = ( params.epoch_build ==~ /default/ ? "" : "-e $params.epoch_build")
        args = ( params.deconv_args.cell2location.build ? params.deconv_args.cell2location.build : "" )
        cuda_device = ( params.gpu ? params.cuda_device : "cpu" )
        println ("Building cell2location model with ${ (params.gpu) ? "GPU" : "CPU" }...")
        """
        source activate cell2loc_env
        python $params.rootdir/subworkflows/deconvolution/cell2location/build_model.py \
            $sc_input $cuda_device -a $params.annot $sample_id_arg $epochs $args -o \$PWD 
        """

}

process fitCell2locationModel {
    tag "c2l_$output_suffix"
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    container 'csangara/sp_cell2location:latest'
    echo true

    input:
        tuple path (sp_input), path (sp_input_rds)
        path (model)
    output:
        tuple val('cell2location'), path("$output"), path (sp_input_rds)
    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_cell2location_${output_suffix}${params.runID_props}.preformat"
        epochs = ( params.epoch_fit ==~ /default/ ? "" : "-e $params.epoch_fit")
        args = ( params.deconv_args.cell2location.fit ? params.deconv_args.cell2location.fit : "" )
        cuda_device = ( params.gpu ? params.cuda_device : "cpu" )
        println ("Fitting cell2location model from file ${model} with ${ (params.gpu) ? "GPU" : "CPU" }...")
        println ("Arguments: ${args}")
        """
        source activate cell2loc_env
        python $params.rootdir/subworkflows/deconvolution/cell2location/fit_model.py \
            $sp_input $model $cuda_device $epochs $args -o \$PWD 
        mv proportions.tsv $output
        """
}
