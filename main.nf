nextflow.enable.dsl=2

include { runMethods } from './scripts/deconvolution/run_methods'
include { computeMetrics } from './scripts/evaluation/evaluate_methods'

workflow {
    main:
        sc_input_ch = Channel.fromPath(params.sc_input)
        sp_input_ch = Channel.fromPath(params.sp_input)

        runMethods(sc_input_ch, sp_input_ch)
        computeMetrics(runMethods.out)
        // propfiles_ch = Channel.fromPath("${params.outdir}/proportions_*${params.output_suffix}")
        // propfiles_ch.view()
        // 

}
