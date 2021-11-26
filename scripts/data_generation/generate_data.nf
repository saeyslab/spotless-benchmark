nextflow.enable.dsl=2
// TODO: ADD REPLICATES AND DATASET TYPES AS INPUT
process generate_synthetic_data {
    echo true
    publishDir "/mnt/d/spade-benchmark/test_generate_data/", mode: 'copy'
    input:
        path (sc_input)
        val (x)

    output:
        path ("${sc_input}_${x}.rds")

    script:
        """
        echo "Received $sc_input"
        echo "hello $x" > ${sc_input}_${x}.rds
        """
}

workflow generateSyntheticData {
    take:
        sc_input
        x
    main:
        generate_synthetic_data(sc_input, x)
    emit:
        generate_synthetic_data.out
}

workflow {
    generateSyntheticData(params.sc_input, "a")
    generateSyntheticData.out.view()
}