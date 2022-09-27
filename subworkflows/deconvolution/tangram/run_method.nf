process getCellComposition {
    tag 'getcellcounts_tangram'
    label 'trivial'
    container 'csangara/seuratdisk:latest'

    input:
        tuple path (sp_input), path(sp_input_rds)

    output:
        tuple path ("composition.csv"), path(sp_input), path(sp_input_rds)
    
    script:
        """
        Rscript $params.rootdir/subworkflows/deconvolution/tangram/getCellComposition.R \
         --sp_input_h5ad $sp_input --sp_input_rds $sp_input_rds
        """

}

process runTangram {
    tag "tangram_${output_suffix}"
    label "retry"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    container 'csangara/sp_tangram:latest'
    echo true

    input:
        path (sc_input)
        tuple path (sp_input), path(sp_input_rds)
        path cell_counts_file // optional file

    output:
        tuple val('tangram'), path("$output"), path (sp_input_rds)

    script:
        args = ( params.deconv_args.tangram ? params.deconv_args.tangram : "" )
        cell_count_arg = cell_counts_file ? "-c $cell_counts_file" : ""
        cuda_device = ( params.gpu ? params.cuda_device : "cpu" )
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_tangram_${output_suffix}${params.runID_props}.preformat"
        println ("Running Tangram with ${ (params.gpu) ? "GPU" : "CPU" }...")
        """
        source activate tangram-env

        python $params.rootdir/subworkflows/deconvolution/tangram/script_nf.py \
            $sc_input $sp_input $cuda_device -a $params.annot -o $output \
            $cell_count_arg $args
        """

}
