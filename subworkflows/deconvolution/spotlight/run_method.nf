process runSpotlight {
    tag "spotlight_$output_suffix"
    label "longer_time"
    label "retry"
    container 'csangara/sp_spotlight:latest'
    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"
    input:
        tuple path (sc_input), path (sp_input)

    output:
        tuple val('spotlight'), path("$output"), path (sp_input)

    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_spotlight_${output_suffix}${params.runID_props}"
        args = (params.deconv_args.spotlight ? params.deconv_args.spotlight : "")
        """
        Rscript $params.rootdir/subworkflows/deconvolution/spotlight/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output $args
        """

}