nextflow.enable.dsl=2

process runMusic {
    tag "$params.output_suffix"
    container 'csangara/spade_music:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        tuple val('music'), file("$output")

    script:
        output = "${params.output}_music${params.output_suffix}"

        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/music/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output --sampleID $params.sampleID
        """

}

process runSpotlight {
    tag "$params.output_suffix"
    container 'csangara/spade_spotlight:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        tuple val('spotlight'), file("$output")

    script:
        output = "${params.output}_spotlight${params.output_suffix}"

        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/spotlight/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output
        """

}

process runRCTD {
    tag "$params.output_suffix"
    container 'csangara/spade_rctd:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        tuple val('rctd'), file("$output")

    script:
        output = "${params.output}_rctd${params.output_suffix}"
        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/rctd/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output
        """

}

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

workflow runMethods {
    // String matching to check which method to run
    main:
        // In input is all, run all methods
        all_methods = "music,rctd,spotlight"
        methods = ( params.methods ==~ /all/ ? all_methods : params.methods )
        if ( methods =~ /music/ ){ runMusic() }
        if ( methods =~ /rctd/ ){ runRCTD() }
        if ( methods =~ /spotlight/ ){ runSpotlight() }

    emit:
        runMusic.out.mix(runRCTD.out, runSpotlight.out)
        
}

