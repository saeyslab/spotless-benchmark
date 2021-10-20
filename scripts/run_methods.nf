nextflow.enable.dsl=2
// DEFAULT PARAMETERS
params.outdir = "/mnt/d/spade-benchmark/deconv_proportions"
params.sc_input = "/mnt/d/spade-benchmark/unit-test/test_sc_data.rds"
params.sp_input = "/mnt/d/spade-benchmark/unit-test/test_sp_data.rds"
params.annot = "subclass"
params.output = "proportions"
params.sampleID = "none"

process runMusic {
    container 'csangara/spade_music:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        file "$output"

    script:
        output = "${params.output}_music"

        """
        Rscript /mnt/d/spade-benchmark/scripts/deconvolution/music/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output --sampleID $params.sampleID
        """

}

process runSpotlight {
    container 'csangara/spade_spotlight:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        file "$output"

    script:
        output = "${params.output}_spotlight"

        """
        Rscript /mnt/d/spade-benchmark/scripts/deconvolution/spotlight/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output
        """

}

process runRCTD {
    container 'csangara/spade_rctd:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        file "$output"

    script:
        output = "${params.output}_rctd"
        """
        Rscript /mnt/d/spade-benchmark/scripts/deconvolution/rctd/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output
        """

}

params.methods = "music,rctd,spotlight"
all_methods = "music,rctd,spotlight"
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
    
}

process computeMetrics {
    container 'csangara/spade_eval:latest'
    script:
        """
        Rscript ../../evaluation/metrics.R
        """
}

//process aggregateMetrics {

//}

workflow {
    main:
        runMethods()
        // propfiles_ch = Channel.fromPath("${params.outdir}/proportions_*")
        // computeMetrics()

}