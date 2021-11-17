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
        epochs = ( params.epoch_build ==~ /default/ ? "" : "-sce $params.epoch_build")

        """
        echo "Received $sc_input, now building stereoscope model..."
        source activate stereoscope
        export LD_LIBRARY_PATH=/opt/conda/envs/stereoscope/lib
        stereoscope run --sc_cnt $sc_input --label_colname $params.annot \
        -n 5000 $epochs -o \$PWD
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
        epochs = ( params.epoch_fit ==~ /default/ ? "" : "-ste $params.epoch_fit")

        """
        echo "Model files $r_file and $logits_file, fitting stereoscope model..."
        source activate stereoscope
        export LD_LIBRARY_PATH=/opt/conda/envs/stereoscope/lib
        stereoscope run --sc_fit $r_file $logits_file \
        --st_cnt $sp_input -n 5000 $epochs -o \$PWD
        mv $sp_file_basename/W*.tsv $output
        """
}


workflow runStereoscope {
    take:
        sc_input
        sp_input
    main:
        /*
        def modeldir = new File("$params.modeldir/stereoscope${params.output_suffix}")

        if (!modeldir.exists()){
            buildStereoscopeModel(convertRDStoH5AD(params.sc_input))
        }
        */

        buildStereoscopeModel(sc_input)
        fitStereoscopeModel(sp_input,
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
        epochs = ( params.epoch_build ==~ /default/ ? "" : "-e $params.epoch_build")
        """
        echo "Building cell2location model..."
        source activate cell2loc_env
        export LD_LIBRARY_PATH=/opt/conda/envs/cell2loc_env/lib
        python $params.rootdir/spade-benchmark/scripts/deconvolution/cell2location/build_model.py \
            $sc_input $params.cuda_device -a $params.annot $sample_id_arg $epochs -o \$PWD 
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
        epochs = ( params.epoch_fit ==~ /default/ ? "" : "-e $params.epoch_fit")

        """
        echo "Model file $model, fitting cell2location model..."
        source activate cell2loc_env
        export LD_LIBRARY_PATH=/opt/conda/envs/cell2loc_env/lib

        python $params.rootdir/spade-benchmark/scripts/deconvolution/cell2location/fit_model.py \
            $sp_input $model $params.cuda_device $epochs -o \$PWD 
        mv proportions.tsv $output
        
        """
}

workflow runCell2location {
    take:
        sc_input
        sp_input
    main:
        buildCell2locationModel(sc_input)
        fitCell2locationModel(sp_input,
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
        output_ch = Channel.empty() // collect output channels

        // R methods
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
        // Python methods
        // First check if there are python methods in the input params
        python_methods = ["stereoscope","cell2location"]
        methods_list = Arrays.asList(methods.split(','))
        if ( !methods_list.disjoint(python_methods) ){

            // Separate pair of inputs into their own channels, then convert
            input = pair_input_ch.multiMap { sc_file, sp_file ->
                                sc_input: sc_file
                                sp_input: sp_file}

            sc_input_h5ad = convert_sc(input.sc_input, params.sc_type)
            sp_input_h5ad = convert_sp(input.sp_input, params.sp_type)

            if ( methods =~ /stereoscope/ ) {
                runStereoscope(sc_input_h5ad, sp_input_h5ad)
                output_ch = output_ch.mix(runStereoscope.out)
            }

            if ( methods =~ /cell2location/ ) {
                runCell2location(sc_input_h5ad, sp_input_h5ad)
                output_ch = output_ch.mix(runCell2location.out)
            }
        }

    emit:
        output_ch
}

