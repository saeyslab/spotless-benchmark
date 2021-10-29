nextflow.enable.dsl=2
include { convertRDStoH5AD as convert_sc ; convertRDStoH5AD as convert_sp } from '../helper_processes'
include { formatTSVFile } from '../helper_processes'

process runMusic {
    tag "$params.output_suffix"
    container 'csangara/spade_music:latest'
    publishDir params.outdir, mode: 'copy'

    output:
        tuple val('music'), path("$output")

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
        tuple val('spotlight'), path("$output")

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
        tuple val('rctd'), path("$output")

    script:
        output = "${params.output}_rctd${params.output_suffix}"
        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/rctd/script_nf.R \
            --sc_input $params.sc_input --sp_input $params.sp_input \
            --annot $params.annot --output $output
        """

}

process buildStereoscopeModel {
    container 'csangara/spade_stereoscope:latest'
    echo true

    input:
        path (sc_input)
    output:
        path "R*.tsv", emit: r_file
        path "logits*.tsv", emit: logits_file

    script:

        """
        echo "Received $sc_input, now building stereoscope model..."
        source activate stereoscope
        export LD_LIBRARY_PATH=/opt/conda/envs/stereoscope/lib
        stereoscope run --sc_cnt $sc_input --label_colname $params.annot \
        -n 1000 -o \$PWD -sce $params.epoch
        """
}
process fitStereoscopeModel {
    container 'csangara/spade_stereoscope:latest'
    echo true
    input:
        path (sp_input)
        path (r_file)
        path (logits_file)
    output:
        tuple val('stereoscope'), path("$output")
    script:
        output = "${params.output}_stereoscope${params.output_suffix}.preformat"
        sp_file_basename = file(sp_input).getSimpleName()

        """
        echo "Model files $r_file and $logits_file, fitting stereoscope model..."
        source activate stereoscope
        export LD_LIBRARY_PATH=/opt/conda/envs/stereoscope/lib
        stereoscope run --sc_fit $r_file $logits_file \
        --st_cnt $sp_input -o \$PWD -n 1000 -ste $params.epoch
        mv $sp_file_basename/W*.tsv $output
        """

}
workflow runStereoscope {
    main:
        /*
        def modeldir = new File("$params.modeldir/stereoscope${params.output_suffix}")

        if (!modeldir.exists()){
            buildStereoscopeModel(convertRDStoH5AD(params.sc_input))
        }
        */

        buildStereoscopeModel(convert_sc(params.sc_input, "seurat"))
        fitStereoscopeModel(convert_sp(params.sp_input, "synthvisium"),
                            buildStereoscopeModel.out.r_file,
                            buildStereoscopeModel.out.logits_file)
        formatTSVFile(fitStereoscopeModel.out)
    emit:
        formatTSVFile.out
    
}

workflow runMethods {
    // String matching to check which method to run
    main:
        // In input is all, run all methods
        all_methods = "music,rctd,spotlight,stereoscope"
        methods = ( params.methods ==~ /all/ ? all_methods : params.methods )
        if ( methods =~ /music/ ){ runMusic() }
        if ( methods =~ /rctd/ ){ runRCTD() }
        if ( methods =~ /spotlight/ ){ runSpotlight() }
        if ( methods =~ /stereoscope/ ) { runStereoscope() }

    emit:
        runMusic.out.mix(runRCTD.out, runSpotlight.out, runStereoscope.out)
        
}

