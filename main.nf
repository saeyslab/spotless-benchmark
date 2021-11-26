nextflow.enable.dsl=2

include { runMethods } from './scripts/deconvolution/run_methods'
include { computeMetrics } from './scripts/evaluation/evaluate_methods'
include { generateSyntheticData } from './scripts/data_generation/generate_data'

bronze_standards = ((1..7)*8).sort().withIndex().collect{ it, i -> "bronze_standard_$it-${ -> i%8+1}".toString() }
all_standards = ["gold_standard_1", "gold_standard_2"] + bronze_standards
all_modes = ["run_standard", "run_dataset", "generate_and_run"]

workflow {
    main:
        if (!all_modes.contains(params.mode)){
            throw new Exception("Error: please enter --mode as 'run_standard' or 'run_dataset'")
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
        // GENERATE SYNTHVISIUM DATA
        else if (params.mode ==~ /generate_and_run/) {
            println("Generating synthvisium data from ${params.sc_input}...")
            // generateSyntheticData(params.sc_input, "a")
            println("Synthetic data copied to ${params.sp_input}")
        }

        // Print inputs (the timing isn't right with with view(), so do this instead)
        // Although view() has the advantage that it gives the absolute path (sc_input_ch.view())
        println("Single-cell reference dataset:")
        println(params.mode ==~ /run_standard/ ? file(params.sc_input)[0]: file(params.sc_input))

        if (params.mode ==~ /generate_and_run/ ){
            // Print the synthvisium files that will be generated
            println("\nSpatial dataset(s) that will be generated:")
            (1..params.synvis_reps).each{
                println "${file(params.sc_input).getSimpleName()}_${params.synvis_type}_rep${it}.rds"
            }
        } // Without glob pattern, there will be multiple files
        else {
            println("\nSpatial dataset(s):") 
            params.sp_input =~ /\*/ ? file(params.sp_input).each{println "$it"} : println (file(params.sp_input))
        }

        sc_input_ch = Channel.fromPath(params.sc_input) // Can only have 1 file
        // Can have one or more files
        // sp_input_ch = (params.mode ==~ /generate_and_run/ ? generateSyntheticData.out : Channel.fromPath(params.sp_input))

        //runMethods(sc_input_ch, sp_input_ch)

        // Files are matched by parsing the names of proportion files
        File sp_path = new File(params.sp_input)
        //computeMetrics(runMethods.out, sp_path.getParent())
        println("End of workflow")
}
