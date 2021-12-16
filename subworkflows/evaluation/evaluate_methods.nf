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