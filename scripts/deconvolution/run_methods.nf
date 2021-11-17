nextflow.enable.dsl=2
include { convertRDStoH5AD as convert_sc ; convertRDStoH5AD as convert_sp } from '../helper_processes'
include { formatTSVFile } from '../helper_processes'

process runMusic {
    tag "$params.output_suffix"
    container 'csangara/spade_music:latest'
    publishDir params.outdir, mode: 'copy'

    input:
        tuple path (sc_input), path (sp_input)
    output:
        tuple val('music'), path("$output")
    
    script:
        output = "${params.output}_music${params.output_suffix}"

        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/music/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output --sampleID $params.sampleID
        """

}

process runSpotlight {
    tag "$params.output_suffix"
    container 'csangara/spade_spotlight:latest'
    publishDir params.outdir, mode: 'copy'

    input:
        tuple path (sc_input), path (sp_input)

    output:
        tuple val('spotlight'), path("$output"), emit: props_file, optional: true

    script:
        output = "${params.output}_spotlight${params.output_suffix}"

        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/spotlight/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output
        """

}

process runRCTD {
    tag "$params.output_suffix"
    container 'csangara/spade_rctd:latest'
    publishDir params.outdir, mode: 'copy'

    input:
        tuple path (sc_input), path (sp_input)

    output:
        tuple val('rctd'), path("$output")

    script:
        output = "${params.output}_rctd${params.output_suffix}"
        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/rctd/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
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
        -n 1000 -o \$PWD -sce $params.epoch_build
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
        --st_cnt $sp_input -o \$PWD -n 1000 -ste $params.epoch_fit
        mv $sp_file_basename/W*.tsv $output
        """
}


workflow runStereoscope {
    take:
        pair_input_ch
    main:
        /*
        def modeldir = new File("$params.modeldir/stereoscope${params.output_suffix}")

        if (!modeldir.exists()){
            buildStereoscopeModel(convertRDStoH5AD(params.sc_input))
        }
        
        buildStereoscopeModel(convert_sc(params.sc_input, "seurat"))
        fitStereoscopeModel(convert_sp(params.sp_input, "synthvisium"),
                            buildStereoscopeModel.out.r_file,
                            buildStereoscopeModel.out.logits_file)
        formatTSVFile(fitStereoscopeModel.out)
        */
        pair_input_ch.multiMap { sc_file, sp_file ->
                                sc_input: sc_file
                                sp_input: sp_file
                                }.set{ input }

        buildStereoscopeModel(input.sc_input)
        fitStereoscopeModel(input.sp_input,
                            buildStereoscopeModel.out.r_file,
                            buildStereoscopeModel.out.logits_file)
        formatTSVFile(fitStereoscopeModel.out)

    emit:
        formatTSVFile.out
    
}

process buildCell2locationModel {
    container 'csangara/spade_cell2location:latest'
    echo true

    input:
        path (sc_input)
    output:
        path "sc.h5ad", emit: model

    script:
        sample_id_arg = ( params.sampleID ==~ /none/ ? "" : "-s $params.sampleID" )
        """
        echo "Building cell2location model..."
        source activate cell2loc_env
        export LD_LIBRARY_PATH=/opt/conda/envs/cell2loc_env/lib
        python $params.rootdir/spade-benchmark/scripts/deconvolution/cell2location/build_model.py \
            $sc_input $params.cuda_device -a $params.annot $sample_id_arg -o \$PWD -e $params.epoch_build
        """

}

process fitCell2locationModel {
    container 'csangara/spade_cell2location:latest'
    echo true

    input:
        path (sp_input)
        path (model)
    output:
        tuple val('cell2location'), path("$output")
    script:
        output = "${params.output}_cell2location${params.output_suffix}.preformat"

        """
        echo "Model file $model, fitting cell2location model..."
        source activate cell2loc_env
        export LD_LIBRARY_PATH=/opt/conda/envs/cell2loc_env/lib

        python $params.rootdir/spade-benchmark/scripts/deconvolution/cell2location/fit_model.py \
            $sp_input $model $params.cuda_device -o \$PWD -e $params.epoch_fit
        mv proportions.tsv $output
        
        """
}

workflow runCell2location {
    take:
        pair_input_ch
    main:
        pair_input_ch.multiMap { sc_file, sp_file ->
                                sc_input: sc_file
                                sp_input: sp_file
                                }.set{ input }

        buildCell2locationModel(input.sc_input)
        fitCell2locationModel(input.sp_input,
                            buildCell2locationModel.out.model)
        formatTSVFile(fitCell2locationModel.out)
    emit:
        formatTSVFile.out
    
}

workflow runMethods {
    take:
        pair_input_ch

    main:
        // String matching to check which method to run
        all_methods = "music,rctd,spotlight,stereoscope,cell2location"
        methods = ( params.methods ==~ /all/ ? all_methods : params.methods )
        output_ch = Channel.empty()

        if ( methods =~ /music/ ){
            runMusic(pair_input_ch)
            output_ch = output_ch.mix(runMusic.out)
        }

        if ( methods =~ /rctd/ ){
            runRCTD(pair_input_ch)
            output_ch = output_ch.mix(runRCTD.out)
        }
        
        if ( methods =~ /spotlight/ ){
            runSpotlight(pair_input_ch)
            output_ch = output_ch.mix(runSpotlight.out)
        }

        if ( methods =~ /stereoscope/ ) {
            runStereoscope(pair_input_ch)
            output_ch = output_ch.mix(runStereoscope.out)
        }

        if ( methods =~ /cell2location/ ) {
            runCell2location(pair_input_ch)
            output_ch = output_ch.mix(runCell2location.out)
        }

    emit:
        output_ch
}

