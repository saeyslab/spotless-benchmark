nextflow.enable.dsl=2

process computeMetrics {
    tag "${method_name}_${output_suffix}"
    container 'csangara/sp_eval:latest'
    publishDir { "${params.outdir.metrics}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy'
    echo true

    input:
        tuple val (method_name), path (props_file), path (sp_input)

    output:
        tuple val (method_name), path ("$metrics_file")

    script:
        output_suffix = file(sp_input).getSimpleName()
        metrics_file = "metrics_${method_name}_${output_suffix}${params.runID_metrics}"
        annot_remap = ( params.remap_annot ? "--remap $params.remap_annot" : "")

        """
        Rscript $params.rootdir/subworkflows/evaluation/metrics.R \
        --props_file $props_file --sp_input $sp_input \
        --output $metrics_file $annot_remap
        echo $metrics_file
        """
}

// In case metrics.R is updated, just rerun metrics computation
// Give sp_input and methods as parameters
workflow {
    sp_input_ch = Channel.fromPath(params.sp_input)
    
    all_methods = "music,rctd,spatialdwls,spotlight,stereoscope,cell2location,destvi,dstg,nnls,seurat,tangram,stride"
    methods = ( params.methods ==~ /all/ ? all_methods : params.methods )
    methods_list = Channel.from(Arrays.asList(methods.split(',')))

    // Get directory of the proportions file
    input_ch = methods_list.combine(sp_input_ch)
                .multiMap { method, sp_input ->
                    inputs: tuple (method, \
                    "${params.outdir.props}/${file(sp_input).getSimpleName().replaceFirst(/_[a-z]{3}[0-9]+/, "")}/proportions_${method}_${file(sp_input).getSimpleName()}${params.runID_props}",\
                    sp_input)
                }
                
    computeMetrics(input_ch.inputs)

}