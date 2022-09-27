process runSpatialDWLS {
    tag "spatialdwls_$output_suffix"
    label "retry"
    container 'csangara/sp_spatialdwls:latest'
    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"
    input:
        tuple path (sc_input), path (sp_input)

    output:
        tuple val('spatialdwls'), path("$output"), path (sp_input)

    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_spatialdwls_${output_suffix}${params.runID_props}"
        args = (params.deconv_args.spatialdwls ? params.deconv_args.spatialdwls : "")
        """
        Rscript $params.rootdir/subworkflows/deconvolution/spatialdwls/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --rootdir $params.rootdir \
            --annot $params.annot --output $output $args
        """

}
