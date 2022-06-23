process runSeurat {
    tag "seurat_$output_suffix"
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2
    container 'csangara/seurat:latest'
    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"

    input:
        tuple path (sc_input), path (sp_input)
    output:
        tuple val('seurat'), path("$output"), path (sp_input)
    
    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_seurat_${output_suffix}${params.runID_props}"
        args = (params.deconv_args.seurat ? params.deconv_args.seurat : "")

        """
        Rscript $params.rootdir/subworkflows/deconvolution/seurat/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output $args
        """

}
