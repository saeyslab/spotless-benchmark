process runRCTD {
    tag "rctd_$output_suffix"
    label "retry"
    container 'csangara/sp_rctd:latest'
    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"
    input:
        tuple path (sc_input), path (sp_input)

    output:
        tuple val('rctd'), path("$output"), path (sp_input)

    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_rctd_${output_suffix}${params.runID_props}"
        num_cores = task.cpus
        args = (params.deconv_args.rctd ? params.deconv_args.rctd : "")

        """
        Rscript $params.rootdir/subworkflows/deconvolution/rctd/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output --num_cores $num_cores $args
        """

}
