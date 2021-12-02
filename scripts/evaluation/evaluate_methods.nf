process computeMetrics {
    container 'csangara/spade_eval:latest'
    publishDir params.outdir.metrics, mode: 'copy'
    tag "$method_name"
    echo true
    input:
        tuple val (method_name), path (props_file), path (sp_input)
    output:
        tuple val (method_name), path ("$metrics_file")
    script:
        // Get corresponding spatial object based on the proportions file name
        sp_filename = file(props_file).getSimpleName().replace("proportions_${method_name}_", "")
        metrics_file = "metrics_${method_name}_${sp_filename}"

        """
        Rscript $params.rootdir/scripts/evaluation/metrics.R \
        --props_file $props_file --sp_input $sp_input --sp_type $params.sp_type \
        --output $metrics_file
        echo $metrics_file
        """
}