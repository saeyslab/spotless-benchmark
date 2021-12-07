nextflow.enable.dsl=2

// Initialization: possible dataset types
synvis_types_map = [aud: "artificial_uniform_distinct", add: "artificial_diverse_distinct",
                    auo: "artificial_uniform_overlap", ado: "artificial_diverse_overlap",
                    adcd: "artificial_dominant_celltype_diverse", apdcd: "artificial_partially_dominant_celltype_diverse",
                    adrcd: "artificial_dominant_rare_celltype_diverse", arrcd: "artificial_regional_rare_celltype_diverse",
                    prior: "prior_from_data"]
synvis_types_fullnames = synvis_types_map.collect{ it.value }
synvis_types_flat = synvis_types_map.collect{[it.key, it.value]}.flatten() // All key and values in a list

process generate_synthetic_data {
    tag "$output"
    container 'csangara/synthvisium:latest'
    publishDir params.outdir.synvis, mode: 'copy'
    input:
        path (sc_input)
        val (dataset_type)
        val (args) // remaining arguments
        each (rep) // run this $rep times

    output:
        path ("$output")

    script:
        output = "${file(sc_input).getSimpleName()}_${dataset_type}_rep${rep}.rds"
        """
        Rscript $params.rootdir/scripts/data_generation/generate_synthetic_data.R \
                --sc_input $sc_input --dataset_type $dataset_type --rep $rep $args
        """
}

workflow generateSyntheticData {
    main:
        // Filter out invalid input, write out full name of abbreviations, filter out duplicates
        synvis_type_input = params.synvis.type.split(',').findAll{ synvis_types_flat.contains(it) }
                                            .collect{ ( synvis_types_fullnames.contains(it) ? it : synvis_types_map[it]) }
                                            .unique()
        // Extra arguments: remove dataset type and number of replicates, then turn to string
        synvis_args_input = params.synvis.findAll{ it.key != "type" && it.key != "reps" && it.key != "sc_input"}
                                  .collect{ "--$it.key $it.value" }.join(" ")

        println("Single-cell reference: $params.synvis.sc_input")
        println("Dataset types to be generated: ${synvis_type_input.join(", ")}")
        println("Number of replicates per dataset type: $params.synvis.reps")
        println("Arguments: ${ (synvis_args_input) ? synvis_args_input: "None (default)" }")

        // Duplicates sc_input by converting it to a list then turning it back into a dataflow
        sc_input_ch = (Channel.fromPath(params.synvis.sc_input).toList()*synvis_type_input.size).flatMap()
        
        generate_synthetic_data(sc_input_ch, Channel.from(synvis_type_input),
                               synvis_args_input, 1..params.synvis.reps.toInteger())
    emit:
        generate_synthetic_data.out
}

workflow {
    generateSyntheticData()
}