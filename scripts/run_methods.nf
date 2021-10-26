// DEFAULT PARAMETERS
params.outdir = "/mnt/d/spade-benchmark/deconv_proportions"
params.sc_input = "/mnt/d/spade-benchmark/unit-test/test_sc_data.rds"
params.sp_input = "/mnt/d/spade-benchmark/unit-test/test_sp_data.rds"
params.sp_type = "synthvisium"
params.annot = "subclass"
params.output = "proportions"
params.output_suffix = ""
params.sampleID = "none"
params.methods = "music,rctd,spotlight"
all_methods = "music,rctd,spotlight"

process runMusic {
    tag "$params.output_suffix"
    container 'csangara/spade_music:latest'
    publishDir params.outdir, mode: 'copy'

    when:
        params.methods =~ /music/ || params.methods ==~ /all/
    output:
        tuple val('music'), file("$output") into music_ch

    script:
        output = "${params.output}_music${params.output_suffix}"

        """
        Rscript /mnt/d/spade-benchmark/scripts/deconvolution/music/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output --sampleID $params.sampleID
        """

}

process runSpotlight {
    tag "$params.output_suffix"
    container 'csangara/spade_spotlight:latest'
    publishDir params.outdir, mode: 'copy'

    when:
        params.methods =~ /spotlight/ || params.methods ==~ /all/
    output:
        tuple val('spotlight'), file("$output") into spotlight_ch

    script:
        output = "${params.output}_spotlight${params.output_suffix}"

        """
        Rscript /mnt/d/spade-benchmark/scripts/deconvolution/spotlight/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output
        """

}

process runRCTD {
    tag "$params.output_suffix"
    container 'csangara/spade_rctd:latest'
    publishDir params.outdir, mode: 'copy'

    when:
        params.methods =~ /rctd/ || params.methods ==~ /all/
    output:
        tuple val('rctd'), file("$output") into rctd_ch

    script:
        output = "${params.output}_rctd${params.output_suffix}"
        """
        Rscript /mnt/d/spade-benchmark/scripts/deconvolution/rctd/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output
        """

}

process computeMetrics {
    container 'csangara/spade_eval:latest'
    echo true
    input:
        tuple val(method_name), file(props_file) from music_ch.mix(
                                                    spotlight_ch,
                                                    rctd_ch)
    output:
        tuple val(method_name), file("$metrics_file") into metrics_ch
    script:
        metrics_file = "metrics_${method_name}${params.output_suffix}"
        """
        echo $method_name
        echo $props_file
        Rscript /mnt/d/spade-benchmark/scripts/evaluation/metrics.R \
        --props_file $props_file --sp_input $params.sp_input --sp_type $params.sp_type \
        --output $metrics_file

        echo $metrics_file
        """
}

// process aggregateMetrics {

// }

/*
process convertRDStoH5AD {
    container 'csangara/seuratdisk:latest'

    script:
        """
        Rscript /mnt/d/spade-benchmark/scripts/deconvolution/convertRDStoH5AD.R \
        --input_path $params.sc_input
        """
}
*/


/*
workflow runMethods {
    // String matching to check which method to run
    main:
        // In input is all, run all methods
        methods = ( params.methods ==~ /all/ ? all_methods : params.methods )
        if( methods =~ /music/ ){
            runMusic()
        }

        if ( methods =~ /rctd/ ){
            runRCTD()
        }
        
        if ( methods =~ /spotlight/ ){
            runSpotlight()
        }
    emit:
        runMusic.out
        runRCTD.out
        runSpotlight.out
        
}

workflow {
    main:
        runMethods()
        Channel.from(runMethods.out).collect().view()
        // convertRDStoH5AD()
        // propfiles_ch = Channel.fromPath("${params.outdir}/proportions_*${params.output_suffix}")
        // propfiles_ch.view()
        // computeMetrics()

}
*/