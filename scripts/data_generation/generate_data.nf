nextflow.enable.dsl=2

params.synvis = [type: "adrcd,blah,artificial_uniform_distinct,aud", reps: 1, clust_var: params.annot]

// n_regions: 5,   clust_var: params.annot, region_var: null, dataset_id: "1",
//                n_spots_min: 50, n_spots_max: 500, visium_mean: 20000, visium_sd: 7000]

synvis_types_map = [aud: "artificial_uniform_distinct", add: "artificial_diverse_distinct",
                    auo: "artificial_uniform_overlap", ado: "artificial_diverse_overlap",
                    adcd: "artificial_dominant_celltype_diverse", apdcd: "artificial_partially_dominant_celltype_diverse",
                    adrcd: "artificial_dominant_rare_celltype_diverse", arrcd: "artificial_regional_rare_celltype_diverse",
                    prior: "prior_from_data"]
synvis_types_fullnames = synvis_types_map.collect{ it.value }
synvis_types_flat = synvis_types_map.collect{[it.key, it.value]}.flatten()

process generate_synthetic_data {
    tag "$output"
    container 'csangara/synthvisium:latest'
    publishDir "${params.rootdir}/spade-benchmark/test_generate_data/", mode: 'copy'
    input:
        path (sc_input)
        val (dataset_type)
        val (args)
        each (rep)

    output:
        path ("$output")

    script:
        output = "${file(sc_input).getSimpleName()}_${dataset_type}_rep${rep}.rds"
        """
        Rscript $params.rootdir/spade-benchmark/scripts/data_generation/generate_synthetic_data.R \
                --sc_input $sc_input --dataset_type $dataset_type --rep $rep $args
        """
}

workflow generateSyntheticData {
    take:
        sc_input
    main:
        // Filter out invalid input, write out full name of abbreviations, only get unique
        synvis_type_input = params.synvis.type.split(',').findAll{ synvis_types_flat.contains(it) }
                                            .collect{ ( synvis_types_fullnames.contains(it) ? it : synvis_types_map[it]) }
                                            .unique()
        // Extra arguments
        synvis_args_input = params.synvis.findAll{ it.key != "type" && it.key != "reps"}
                                  .collect{ "--$it.key $it.value" }.join(" ")
        println("Single-cell reference: $sc_input")
        println("Dataset types to be generated: ${synvis_type_input.join(", ")}")
        println("Number of replicates per dataset type: $params.synvis.reps")
        println("Arguments: ${ (synvis_args_input) ? synvis_args_input: "None (default)" }")
                             
        generate_synthetic_data(sc_input, Channel.from(synvis_type_input),
                                synvis_args_input, 1..params.synvis.reps.toInteger())
    //emit:
    //    generate_synthetic_data.out
}

workflow {
    generateSyntheticData(params.sc_input)
}