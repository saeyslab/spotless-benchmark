nextflow.enable.dsl=2

// Deconvolution methods
include { runMusic } from './music/run_method.nf'
include { runRCTD } from './rctd/run_method.nf'
include { runSpotlight } from './spotlight/run_method.nf'
include { runSpatialDWLS } from './spatialdwls/run_method.nf'
include { buildCell2locationModel; fitCell2locationModel} from './cell2location/run_method.nf'
include { buildStereoscopeModel; fitStereoscopeModel } from './stereoscope/run_method.nf'
include { buildDestVIModel; fitDestVIModel } from './destvi/run_method.nf'
include { runDSTG } from './dstg/run_method.nf'
include { runNNLS } from './nnls/run_method.nf'
include { runSeurat } from './seurat/run_method.nf'
include { getCellComposition; runTangram } from './tangram/run_method.nf'
include { buildSTRIDEModel; fitSTRIDEModel; runSTRIDE } from './stride/run_method.nf'

// Helper functions
include { convertBetweenRDSandH5AD as convert_sc ; convertBetweenRDSandH5AD as convert_sp } from '../helper_processes'
include { formatTSVFile as formatStereoscope; formatTSVFile as formatC2L; formatTSVFile as formatDestVI; formatTSVFile as formatDSTG;
          formatTSVFile as formatTangram; formatTSVFile as formatSTRIDE } from '../helper_processes'
include { createDummyFile } from '../helper_processes'

workflow runMethods {
    take:
        sc_input_ch
        sp_input_ch
        sc_input_type
        sp_input_type

    main:
        // String matching to check which method to run
        all_methods = "music,rctd,spatialdwls,spotlight,stereoscope,cell2location,destvi,dstg,nnls,seurat,tangram,stride"
        r_methods = ["music", "rctd", "spatialdwls", "spotlight", "dstg", "nnls", "seurat"]
        python_methods = ["stereoscope","cell2location","destvi","tangram","stride"]

        methods = ( params.methods.toLowerCase() ==~ /all/ ? all_methods : params.methods.toLowerCase() )
        methods_list = Arrays.asList(methods.split(','))
        
        output_ch = Channel.empty() // collect output channels

        // R methods
        if ( !methods_list.disjoint(r_methods) ){
            // If the input file is H5AD, convert to RDS
            sc_input_R = sc_input_type ==~ /rds/ ? sc_input_ch :
                                                   convert_sc(sc_input_ch).flatten().filter( ~/.*rds*/ )
            sp_input_R = sp_input_type ==~ /rds/ ? sp_input_ch :
                                                   convert_sp(sp_input_ch).flatten().filter( ~/.*rds*/ )

            pair_input_ch = sc_input_R.combine(sp_input_R)

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

            if ( methods =~ /spatialdwls/ ){
                runSpatialDWLS(pair_input_ch)
                output_ch = output_ch.mix(runSpatialDWLS.out)
            }

            if ( methods =~ /dstg/ ){
                runDSTG(pair_input_ch)
                formatDSTG(runDSTG.out) 
                output_ch = output_ch.mix(formatDSTG.out)
            }

            if ( methods =~ /nnls/ ){
                runNNLS(pair_input_ch)
                output_ch = output_ch.mix(runNNLS.out)
            }

            if ( methods =~ /seurat/ ){
                runSeurat(pair_input_ch)
                output_ch = output_ch.mix(runSeurat.out)
            }
        }
        // Python methods
        // First check if there are python methods in the input params
        // before performing conversion of data to h5ad
        
        if ( !methods_list.disjoint(python_methods) ){
            // If the input file is RDS, convert to H5AD
            sc_input_conv = sc_input_type ==~ /rds/ ? convert_sc(sc_input_ch).flatten().filter( ~/.*h5ad*/ ) :
                                                      sc_input_ch

            // If the spatial file is H5AD, create dummy file
            sp_input_pair = sp_input_type ==~ /rds/ ? 
                            convert_sp(sp_input_ch) : createDummyFile(sp_input_ch)

            if ( methods =~ /stereoscope/ ) {
                buildStereoscopeModel(sc_input_conv)

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
                buildCell2locationModel(sc_input_conv)

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
                
                buildDestVIModel(sc_input_conv, sp_input_pair.first())
                
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

            if ( methods =~ /tangram/ ){
                // For constrained mode, need to run extra function
                if (params.deconv_args.tangram =~ /-m[ ]*constrained/) {
                    // Get cell counts from spatial file (or segment)
                    getCellComposition(sp_input_pair)

                    // Repeat cell composition and spatial file for each single-cell file
                    getCellComposition.out.combine(sc_input_conv)
                    .multiMap { cell_count_file, sp_file_h5ad, sp_file_rds, sc_file ->
                                sc_input: sc_file
                                sp_input: tuple sp_file_h5ad, sp_file_rds
                                cell_count: cell_count_file }
                    .set{ tangram_combined_ch }
                    
                    runTangram(tangram_combined_ch.sc_input,  tangram_combined_ch.sp_input,
                               tangram_combined_ch.cell_count)
                } else {
                    // Only take single-cell and spatial file
                    sc_input_conv.combine(sp_input_pair)
                    .multiMap { sc_file, sp_file_h5ad, sp_file_rds ->
                                sc_input: sc_file
                                sp_input: tuple sp_file_h5ad, sp_file_rds}
                    .set{ tangram_combined_ch }

                    // Cell count file as empty list
                    runTangram(tangram_combined_ch.sc_input,  tangram_combined_ch.sp_input, [])

                }

                formatTangram(runTangram.out)
                output_ch = output_ch.mix(formatTangram.out)

            }

            if ( methods =~ /stride/ ) {
                
                /*
                buildSTRIDEModel(sc_input_conv)

                // Repeat model output for each spatial file
                buildSTRIDEModel.out.combine(sp_input_pair)
                .multiMap { zip_file, txt_files, sp_file_h5ad, sp_file_rds ->
                            model: tuple zip_file, txt_files
                            sp_input: tuple sp_file_h5ad, sp_file_rds }
                .set{ stride_combined_ch }
                
                fitSTRIDEModel(stride_combined_ch.sp_input,
                               stride_combined_ch.model)
                formatSTRIDE(fitSTRIDEModel.out) 
                */
                
                
                runSTRIDE(sc_input_conv.combine(sp_input_pair))
                formatSTRIDE(runSTRIDE.out)
                
                output_ch = output_ch.mix(formatSTRIDE.out)
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
