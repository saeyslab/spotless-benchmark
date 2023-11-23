nextflow.enable.dsl=2

include { runMethods } from './subworkflows/deconvolution/run_methods'
include { computeMetrics } from './subworkflows/evaluation/evaluate_methods'
include { generateSyntheticData } from './subworkflows/data_generation/generate_data'

silver_standards = ((1..7)*9).sort().withIndex().collect{ it, i -> "silver_standard_$it-${ -> i%9+1}".toString() }
silver_standards_coarse = (1..7).collect{ it -> "silver_standard_$it".toString()}
gold_standards = (1..3).collect{ it -> "gold_standard_$it".toString()}
all_standards = gold_standards + silver_standards + silver_standards_coarse + ["silver_standard_1-10", "silver_standard_1-11"]
all_modes = ["run_standard", "run_dataset", "generate_and_run"]


workflow {
    main:
               
        if (!all_modes.contains(params.mode)){
            throw new Exception("Error: please enter --mode as 'run_standard', 'run_dataset', or 'generate_and_run'")
        }

        // RUN STANDARD PIPELINE
        if (params.mode ==~ /run_standard/ ){
            // Check pipeline name
            if (!all_standards.contains(params.standard)){
                throw new Exception("Error: standard '$params.standard' does not exist.")
            }
            // Check if it was run with config file
            if (params.sc_input =~ /test_sc_data/ && params.sp_input =~ /test_sp_data/){
                throw new Exception("Error: default input detected, please run with '-c data/pipeline.config'")
            }
            println("Running the pipeline on $params.standard...")
        } 
        // RUN ON YOUR OWN DATA
        else if (params.mode ==~ /run_dataset/ ) {
            println("Running the pipeline on the provided data...")
        }
        // GENERATE SYNTHETIC DATA AND RUN
        else if (params.mode ==~ /generate_and_run/) {
            println("Generating synthetic data with synthspot...")
            generateSyntheticData()
            println("Synthetic data will be copied to ${params.outdir.synthspot}")
        }

        // Print inputs (the timing isn't right with with view(), so do this instead)
        // Although view() has the advantage that it gives the absolute path (sc_input_ch.view())
        if (!(params.mode ==~ /generate_and_run/)) {
            if (params.verbose) {
                println("Single-cell reference dataset:")
                println(params.mode ==~ /run_standard/ ? file(params.sc_input)[0]: file(params.sc_input))

                println("\nSpatial dataset(s):") 
                // With glob pattern, there will be multiple files
                params.sp_input =~ /\*/ ? file(params.sp_input).each{println "$it"} : println (file(params.sp_input))
            }
        }

        sc_input_type = params.sc_input.split("\\.").last() =~ /h5/ ? "h5ad" : "rds"
        sp_input_type = params.sp_input.split("\\.").last() =~ /h5/ ? "h5ad" : "rds"
        
        println("The single-cell data is of ${sc_input_type} format and the spatial data is of ${sp_input_type} format.")

        // Print arguments
        println("\nParameters:")
        params.deconv_args.each{ method, argu  -> println "${method}: ${argu}" }

        if (sp_input_type ==~ /h5ad/ && !params.skip_metrics){
            throw new Exception("Metric calculation is not supported for spatial h5ad files as input, please use --skip_metrics if you wish to run the pipeline.")
        }

        // Prepare the input
        sc_input_ch = Channel.fromPath(params.sc_input) // Can only have 1 file
        // Can have one or more files
        sp_input_ch = (params.mode ==~ /generate_and_run/ ? generateSyntheticData.out : Channel.fromPath(params.sp_input))

        runMethods(sc_input_ch, sp_input_ch, sc_input_type, sp_input_type)
        
        if (!params.skip_metrics){
            computeMetrics(runMethods.out)
        }

}
