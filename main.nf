nextflow.enable.dsl=2

include { runMethods } from './scripts/deconvolution/run_methods'
include { computeMetrics } from './scripts/evaluation/evaluate_methods'
include { generate_synthetic_data } from './scripts/data_generation/generate_data'

bronze_standards = ((1..7)*8).sort().withIndex().collect{ it, i -> "bronze_standard_$it-${ -> i%8+1}".toString() }
all_pipelines = ["gold_standard_1", "gold_standard_2"] + bronze_standards

workflow {
    main:
        // Run standard
        if (!( params.pipeline ==~ /none/ )) {
            if (!all_pipelines.contains(params.pipeline)){
                throw new Exception("Error: pipeline not found.")
            }
            println(params.pipeline)
        }
        /*
        if ( params.mode =~ /generate/ ) {
            generate_synthetic_data(params.sc_input)

            if ( params.mode ==~ /generate_and_run_methods/ ) {
                sp_input_ch = Channel.fromPath(params.sp_input)
            }
        }

    
        sc_input_ch = Channel.fromPath(params.sc_input) // Can have 1 file
        sp_input_ch = Channel.fromPath(params.sp_input) // Can have one or more files

        runMethods(sc_input_ch, sp_input_ch)

        // Files are matched by parsing the names of proportion files
        File sp_path = new File(params.sp_input)
        computeMetrics(runMethods.out, sp_path.getParent())
        */
        println("hello world")
}
