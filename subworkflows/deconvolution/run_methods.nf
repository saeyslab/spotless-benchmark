nextflow.enable.dsl=2

// Deconvolution methods
include { runMusic } from './music/run_method.nf'
include { runRCTD } from './rctd/run_method.nf'
include { runSpotlight } from './spotlight/run_method.nf'
include { buildCell2locationModel; fitCell2locationModel} from './cell2location/run_method.nf'
include { buildStereoscopeModel; fitStereoscopeModel } from './stereoscope/run_method.nf'
include { buildDestVIModel; fitDestVIModel } from './destvi/run_method.nf'

// Helper functions
include { convertRDStoH5AD as convert_sc ; convertRDStoH5AD as convert_sp } from '../helper_processes'
include { formatTSVFile as formatStereoscope; formatTSVFile as formatC2L; formatTSVFile as formatDestVI } from '../helper_processes'


workflow runMethods {
    take:
        sc_input_ch
        sp_input_ch

    main:
        // String matching to check which method to run
        all_methods = "music,rctd,spotlight,stereoscope,cell2location,destvi"
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
        python_methods = ["stereoscope","cell2location","destvi"]
        methods_list = Arrays.asList(methods.split(','))
        if ( !methods_list.disjoint(python_methods) ){

            sc_input_pair = convert_sc(sc_input_ch)
            sp_input_pair = convert_sp(sp_input_ch)

            if ( methods =~ /stereoscope/ ) {
                buildStereoscopeModel(sc_input_pair)

                // Repeat model output for each spatial file
                buildStereoscopeModel.out.combine(sp_input_pair)
                .multiMap { r_file, logits_file, sp_file_h5ad, sp_file_rds ->
                            model: tuple r_file, logits_file
                            sp_input: tuple sp_file_h5ad, sp_file_rds }
                .set{ stereo_combined_ch }

                fitStereoscopeModel(stereo_combined_ch.sp_input,
                                    stereo_combined_ch.model)
                formatStereoscope(fitStereoscopeModel.out) 
                output_ch = output_ch.mix(formatStereoscope.out)
            }

            if ( methods =~ /cell2location/ ) {
                buildCell2locationModel(sc_input_pair)

                // Repeat model output for each spatial file
                buildCell2locationModel.out.combine(sp_input_pair)
                .multiMap { model_sc_file, sp_file_h5ad, sp_file_rds ->
                            model: model_sc_file
                            sp_input: tuple sp_file_h5ad, sp_file_rds }
                .set{ c2l_combined_ch }

                fitCell2locationModel(c2l_combined_ch.sp_input,
                                      c2l_combined_ch.model)
                formatC2L(fitCell2locationModel.out)
                output_ch = output_ch.mix(formatC2L.out)
            }

            if ( methods =~ /destvi/ ) {
                buildDestVIModel(sc_input_pair)

                // Repeat model output for each spatial file
                buildDestVIModel.out.combine(sp_input_pair)
                .multiMap { h5ad_file, model_file, sp_file_h5ad, sp_file_rds ->
                            model: tuple h5ad_file, model_file
                            sp_input: tuple sp_file_h5ad, sp_file_rds }
                .set{ destvi_combined_ch }

                fitDestVIModel(destvi_combined_ch.sp_input,
                               destvi_combined_ch.model)
                formatDestVI(fitDestVIModel.out) 
                output_ch = output_ch.mix(formatDestVI.out)
            }
        }

    emit:
        output_ch
}

workflow {

    if (!(params.mode ==~ /run_dataset/)){
        throw new Exception("Error: can only run this with the 'run_dataset' mode")
    }

    // RUN ON YOUR OWN DATA
    println("Running the pipeline on the provided data...")
    
    // Print inputs (the timing isn't right with with view(), so do this instead)
    // Although view() has the advantage that it gives the absolute path (sc_input_ch.view())
    if (params.verbose) {
        println("Single-cell reference dataset:")
        println(file(params.sc_input))

        println("\nSpatial dataset(s):") 
        // With glob pattern, there will be multiple files
        params.sp_input =~ /\*/ ? file(params.sp_input).each{println "$it"} : println (file(params.sp_input))
    }

    sc_input_ch = Channel.fromPath(params.sc_input) // Can only have 1 file
    sp_input_ch = Channel.fromPath(params.sp_input) // Can have one or more files

    runMethods(sc_input_ch, sp_input_ch)
}
