nextflow.enable.dsl=2

include { runMethods } from './scripts/deconvolution/run_methods'
include { computeMetrics } from './scripts/evaluation/evaluate_methods'
include { generateSyntheticData } from './scripts/data_generation/generate_data'

bronze_standards = ((1..7)*8).sort().withIndex().collect{ it, i -> "bronze_standard_$it-${ -> i%8+1}".toString() }
all_pipelines = ["gold_standard_1", "gold_standard_2"] + bronze_standards
all_modes = ["run_pipeline", "generate_and_run", "run_own_dataset"]

workflow {
    main:
        if (!all_modes.contains(params.mode)){
            throw new Exception("Error: please enter --mode as 'run_pipeline', 'generate_and_run' or 'run_own_dataset'")
        }

        // RUN STANDARD PIPELINE
        if (params.mode ==~ /run_pipeline/ ){
            // Check pipeline name
            if (!all_pipelines.contains(params.pipeline)){
                throw new Exception("Error: pipeline '$params.pipeline' does not exist.")
            }
            // Check if it was run with config file
            if (params.sc_input =~ /test_sc_data/ && params.sp_input =~ /test_sp_data/){
                throw new Exception("Error: default input detected, please run with '-c data/pipeline.config'")
            }
            println("Running $params.pipeline pipeline...")
        }
        
        // GENERATE SYNTHVISIUM DATA
        else if (params.mode ==~ /generate_and_run/) {
            println("Generating synthvisium data from ${params.sc_input}...")
            // generateSyntheticData(params.sc_input, params.sp_input)
            println("Synthetic data saved at ${params.sp_input}")
        }

        // RUN ON YOUR OWN DATA
        else {
            println("Running the pipeline on the provided data...")
        }
        
        sc_input_ch = Channel.fromPath(params.sc_input) // Can have 1 file
        sp_input_ch = Channel.fromPath(params.sp_input) // Can have one or more files

        // Print inputs (the timing isn't right with with view(), so do this instead)
        // However, view() has the advantage that it gives the absolute path (sc_input_ch.view())
        println("Single-cell reference dataset:")
        println(params.mode ==~ /run_pipeline/ ? file(params.sc_input)[0]: file(params.sc_input))
        println("\nSpatial dataset(s):")
        params.sp_input =~ /\*/ ? file(params.sp_input).each{println "$it"} : println (file(params.sp_input))

        //runMethods(sc_input_ch, sp_input_ch)

        // Files are matched by parsing the names of proportion files
        File sp_path = new File(params.sp_input)
        //computeMetrics(runMethods.out, sp_path.getParent())
        println("End of workflow")
}
