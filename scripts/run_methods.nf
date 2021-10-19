nextflow.enable.dsl=2
params.outdir = "/mnt/d/spade-benchmark/deconv_proportions"
process runMusic {
    container 'csangara/spade_music:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        file 'proportions_music'

    script:
        """
        Rscript ../../deconvolution/music/script_nf.R
        """

}

process runSpotlight {
    container 'csangara/spade_spotlight:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        file 'proportions_spotlight'

    script:
        """
        Rscript ../../deconvolution/spotlight/script_nf.R
        """

}

process runRCTD {
    container 'csangara/spade_rctd:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        file 'proportions_rctd'

    script:
        """
        Rscript ../../deconvolution/rctd/script_nf.R
        """

}

params.methods = "music,rctd,spotlight"

workflow runMethods {
    // String matching to check which method to run
    main:
        if( params.methods =~ /music/ ){
            runMusic()
        }

        if ( params.methods =~ /rctd/ ){
            runRCTD()
        }
        
        if ( params.methods =~ /spotlight/ ){
            runSpotlight()
        }
    
}

workflow {
    main:
        runMethods()
        propfiles_ch = Channel.fromPath("${params.outdir}/proportions_*")
        propfiles_ch.view()
}