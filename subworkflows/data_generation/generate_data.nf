nextflow.enable.dsl=2
// Helper functions
include { convertBetweenRDSandH5AD as convert_sc } from '../helper_processes'

// Initialization: possible dataset types
synthspot_types_map = [aud: "artificial_uniform_distinct", add: "artificial_diverse_distinct",
                    auo: "artificial_uniform_overlap", ado: "artificial_diverse_overlap",
                    adcd: "artificial_dominant_celltype_diverse", apdcd: "artificial_partially_dominant_celltype_diverse",
                    adrcd: "artificial_dominant_rare_celltype_diverse", arrcd: "artificial_regional_rare_celltype_diverse",
                    prior: "prior_from_data", real: "real", rm: "real_missing_celltypes_visium",
                    arm: "artificial_missing_celltypes_visium", addm: "artificial_diverse_distinct_missing_celltype_sc",
                    adom: "artificial_diverse_overlap_missing_celltype_sc"]
synthspot_types_fullnames = synthspot_types_map.collect{ it.value }
synthspot_types_flat = synthspot_types_map.collect{[it.key, it.value]}.flatten() // All key and values in a list

process generate_synthetic_data {
    tag "$output"
    container 'csangara/synthspot:latest'
    publishDir params.outdir.synthspot, mode: 'copy'
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
        Rscript $params.rootdir/subworkflows/data_generation/generate_synthetic_data.R \
                --sc_input $sc_input --dataset_type $dataset_type --rep $rep $args
        """
}

workflow generateSyntheticData {
    main:
        // Check input type - if h5ad, convert to rds
        sc_input_type = file(params.synthspot.sc_input).getExtension() =~ /h5/ ? "h5ad" : "rds"
        println("The synthetic data is of ${sc_input_type} format.")
        sc_input_ch = Channel.fromPath(params.synthspot.sc_input)
        sc_input_conv = sc_input_type =~ /h5ad/ ? convert_sc(sc_input_ch).flatten().filter( ~/.*rds*/ ) : sc_input_ch
        
        // Filter out invalid input, write out full name of abbreviations, filter out duplicates
        synthspot_type_input = params.synthspot.type.split(',').findAll{ synthspot_types_flat.contains(it) }
                                            .collect{ ( synthspot_types_fullnames.contains(it) ? it : synthspot_types_map[it]) }
                                            .unique()
        // Extra arguments: remove dataset type and number of replicates, then turn to string
        synthspot_args_input = params.synthspot.findAll{ it.key != "type" && it.key != "reps" && it.key != "sc_input"}
                                  .collect{ "--$it.key $it.value" }.join(" ")

        println("Single-cell reference: $params.synthspot.sc_input")
        println("Dataset types to be generated: ${synthspot_type_input.join(", ")}")
        println("Number of replicates per dataset type: $params.synthspot.reps")
        println("Arguments: ${ (synthspot_args_input) ? synthspot_args_input: "None (default)" }")

        // Duplicates sc_input by converting it to a list then turning it back into a dataflow
        sc_input_dup = (sc_input_conv.toList()*synthspot_type_input.size).flatMap()
        
        generate_synthetic_data(sc_input_dup, Channel.from(synthspot_type_input),
                               synthspot_args_input, 1..params.synthspot.reps.toInteger())
    emit:
        generate_synthetic_data.out
}

workflow {
    generateSyntheticData()
}