nextflow.enable.dsl=2
include { convertRDStoH5AD as convert_sc ; convertRDStoH5AD as convert_sp } from '../helper_processes'
include { formatTSVFile as formatStereoscope; formatTSVFile as formatC2L } from '../helper_processes'

process runMusic {
    tag "music_$output_suffix"
    container 'csangara/spade_music:latest'
    publishDir params.outdir.props, mode: 'copy'

    input:
        tuple path (sc_input), path (sp_input)
    output:
        tuple val('music'), path("$output"), path (sp_input)
    
    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_music_${output_suffix}"

        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/music/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output --sampleID $params.sampleID
        """

}

process runSpotlight {
    tag "spotlight_$output_suffix"
    container 'csangara/spade_spotlight:latest'
    publishDir params.outdir.props, mode: 'copy'

    input:
        tuple path (sc_input), path (sp_input)

    output:
        tuple val('spotlight'), path("$output"), path (sp_input)

    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_spotlight_${output_suffix}"
        args = (params.deconv_args.spotlight ? params.deconv_args.spotlight : "")
        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/spotlight/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output $args
        """

}

process runRCTD {
    tag "rctd_$output_suffix"
    container 'csangara/spade_rctd:latest'
    publishDir params.outdir.props, mode: 'copy'

    input:
        tuple path (sc_input), path (sp_input)

    output:
        tuple val('rctd'), path("$output"), path (sp_input)

    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_rctd_${output_suffix}"
        """
        Rscript $params.rootdir/spade-benchmark/scripts/deconvolution/rctd/script_nf.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output
        """

}

process buildStereoscopeModel {
    tag 'stereo_build'
    label "retry"
    container 'csangara/spade_stereoscope:latest'
    echo true

    input:
        path (sc_input)
    output:
        tuple path ("R*.tsv"), path ("logits*.tsv")

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
    tag "stereo_$sp_file_basename"
    label "retry"
    container 'csangara/spade_stereoscope:latest'
    echo true

    input:
        path (sp_input)
        tuple path (r_file), path (logits_file)
    output:
        tuple val('stereoscope'), path("$output"), path (sp_input)
    script:
        sp_file_basename = file(sp_input).getSimpleName()
        output = "proportions_stereoscope_${sp_file_basename}.preformat"
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

process buildCell2locationModel {
    tag 'c2l_build'
    label "retry"
    container 'csangara/spade_cell2location:latest'
    echo true

    input:
        path (sc_input)
    output:
        path "sc.h5ad"

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
    tag "c2l_$output_suffix"
    label "retry"
    container 'csangara/spade_cell2location:latest'
    echo true

    input:
        path (sp_input)
        path (model)
    output:
        tuple val('cell2location'), path("$output"), path (sp_input)
    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_cell2location_${output_suffix}.preformat"
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

workflow runMethods {
    take:
        sc_input_ch
        sp_input_ch

    main:
        // String matching to check which method to run
        all_methods = "music,rctd,spotlight,stereoscope,cell2location"
        methods = ( params.methods ==~ /all/ ? all_methods : params.methods )
        output_ch = Channel.empty() // collect output channels

        // R methods
        pair_input_ch = sc_input_ch.combine(sp_input_ch)
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
        // before performing conversion of data to h5ad
        python_methods = ["stereoscope","cell2location"]
        methods_list = Arrays.asList(methods.split(','))
        if ( !methods_list.disjoint(python_methods) ){

            sc_input_h5ad = convert_sc(sc_input_ch, params.sc_type)
            sp_input_h5ad = convert_sp(sp_input_ch, params.sp_type)

            if ( methods =~ /stereoscope/ ) {
                buildStereoscopeModel(sc_input_h5ad)

                // Repeat model output for each spatial file
                buildStereoscopeModel.out.combine(sp_input_h5ad)
                .multiMap { r_file, logits_file, sp_file ->
                            model: tuple r_file, logits_file
                            sp_input_h5ad: sp_file }
                .set{ stereo_combined_ch }

                fitStereoscopeModel(stereo_combined_ch.sp_input_h5ad,
                                    stereo_combined_ch.model)
                formatStereoscope(fitStereoscopeModel.out) 
                output_ch = output_ch.mix(formatStereoscope.out)
            }

            if ( methods =~ /cell2location/ ) {
                buildCell2locationModel(sc_input_h5ad)

                // Repeat model output for each spatial file
                buildCell2locationModel.out.combine(sp_input_h5ad)
                .multiMap { model_sc_file, sp_file ->
                            model: model_sc_file
                            sp_input_h5ad: sp_file }
                .set{ c2l_combined_ch }

                fitCell2locationModel(c2l_combined_ch.sp_input_h5ad,
                                      c2l_combined_ch.model)
                formatC2L(fitCell2locationModel.out)
                output_ch = output_ch.mix(formatC2L.out)
            }
        }

    emit:
        output_ch
}

