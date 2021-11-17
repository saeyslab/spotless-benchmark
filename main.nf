nextflow.enable.dsl=2

include { runMethods } from './scripts/deconvolution/run_methods'
include { computeMetrics } from './scripts/evaluation/evaluate_methods'

workflow {
    main:
        sc_input_ch = Channel.fromPath(params.sc_input)
        sp_input_ch = Channel.fromPath(params.sp_input)
        pair_input_ch = sc_input_ch.combine(sp_input_ch)

        runMethods(pair_input_ch)
        computeMetrics(runMethods.out)
        // propfiles_ch = Channel.fromPath("${params.outdir}/proportions_*${params.output_suffix}")
        // propfiles_ch.view()
        // 

}
