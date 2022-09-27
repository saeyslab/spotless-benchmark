process runMusic {
    tag "music_$output_suffix"
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2
    container 'csangara/sp_music:latest'
    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"

    input:
        tuple path (sc_input), path (sp_input)
    output:
        tuple val('music'), path("$output"), path (sp_input)
    
    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_music_${output_suffix}${params.runID_props}"
        args = (params.deconv_args.music ? params.deconv_args.music : "")

        """
        Rscript $params.rootdir/subworkflows/deconvolution/music/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output \
            --rootdir $params.rootdir \
            --sampleID $params.sampleID $args
        """

}
