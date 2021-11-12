nextflow.enable.dsl=2

include { runMethods } from './scripts/deconvolution/run_methods'
include { computeMetrics } from './scripts/evaluation/evaluate_methods'


workflow {
    main:
        runMethods()
        computeMetrics(runMethods.out)
        // convertRDStoH5AD()
        // propfiles_ch = Channel.fromPath("${params.outdir}/proportions_*${params.output_suffix}")
        // propfiles_ch.view()
        // 

}