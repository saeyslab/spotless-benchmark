nextflow.enable.dsl=2

process computeMetrics {
    tag "$method_name"
    container 'csangara/spade_eval:latest'
    publishDir { "${params.outdir.metrics}/${output_suffix.replaceFirst(/_rep[0-9]+/, "")}" },
                mode: 'copy'
    echo true

    input:
        tuple val (method_name), path (props_file), path (sp_input)

    output:
        tuple val (method_name), path ("$metrics_file")

    script:
        output_suffix = file(sp_input).getSimpleName()
        metrics_file = "metrics_${method_name}_${output_suffix}"

        """
        Rscript $params.rootdir/subworkflows/evaluation/metrics.R \
        --props_file $props_file --sp_input $sp_input --sp_type $params.sp_type \
        --output $metrics_file
        echo $metrics_file
        """
}

// In case metrics.R is updated, just rerun metrics computation
// Give sp_input and methods as parameters
workflow {
    sp_input_ch = Channel.fromPath(params.sp_input)
    
    all_methods = "music,rctd,spotlight,stereoscope,cell2location"
    methods = ( params.methods ==~ /all/ ? all_methods : params.methods )
    methods_list = Channel.from(Arrays.asList(methods.split(',')))

    // Get directory of the proportions file
    input_ch = methods_list.combine(sp_input_ch)
                .multiMap { method, sp_input ->
                    inputs: tuple (method, \
                    "${params.outdir.props}/${file(sp_input).getSimpleName().replaceFirst(/_rep[0-9]+/, "")}/proportions_${method}_${file(sp_input).getSimpleName()}",\
                    sp_input)
                }
                
    computeMetrics(input_ch.inputs)

}