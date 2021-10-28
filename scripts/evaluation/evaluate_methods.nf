
/*
params.sp_input = "/mnt/d/spade-benchmark/unit-test/test_sp_data.rds"
params.sp_type = "synthvisium"
params.output_suffix = ""
*/
process computeMetrics {
    container 'csangara/spade_eval:latest'
    echo true
    input:
        tuple val(method_name), file(props_file)
    output:
        tuple val(method_name), file("$metrics_file")
    script:
        metrics_file = "metrics_${method_name}${params.output_suffix}"
        """
        echo $method_name
        echo $props_file
        # Rscript $params.rootdir/spade-benchmark/scripts/evaluation/metrics.R \
        # --props_file $props_file --sp_input $params.sp_input --sp_type $params.sp_type \
        # --output $metrics_file

        echo 'hello world' > $metrics_file
        """
}