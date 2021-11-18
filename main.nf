nextflow.enable.dsl=2

include { runMethods } from './scripts/deconvolution/run_methods'
include { computeMetrics } from './scripts/evaluation/evaluate_methods'

workflow {
    main:
        sc_input_ch = Channel.fromPath(params.sc_input) // Can have 1 file
        sp_input_ch = Channel.fromPath(params.sp_input) // Can have one or more files

        runMethods(sc_input_ch, sp_input_ch)

        // Files are matched by parsing the names of proportion files
        File sp_path = new File(params.sp_input)
        computeMetrics(runMethods.out, sp_path.getParent())

}
