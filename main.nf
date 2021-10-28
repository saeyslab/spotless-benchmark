nextflow.enable.dsl=2

include { runMethods } from './scripts/deconvolution/run_methods'
include { computeMetrics } from './scripts/evaluation/evaluate_methods'

params.outdir = "/mnt/d/spade-benchmark/deconv_proportions"
params.sc_input = "/mnt/d/spade-benchmark/unit-test/test_sc_data.rds"
params.sp_input = "/mnt/d/spade-benchmark/unit-test/test_sp_data.rds"
params.sp_type = "synthvisium"
params.annot = "subclass"
params.output = "proportions"
params.output_suffix = ""
params.sampleID = "none"
params.methods = "music,rctd,spotlight"


workflow {
    main:
        runMethods()
        computeMetrics(runMethods.out)
        // convertRDStoH5AD()
        // propfiles_ch = Channel.fromPath("${params.outdir}/proportions_*${params.output_suffix}")
        // propfiles_ch.view()
        // 

}