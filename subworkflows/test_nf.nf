nextflow.enable.dsl=2


process foo {
  echo true
  input:
    tuple path (sc_input), path (sp_input)
  script:
    """
    echo "sc: $sc_input, sp: $sp_input"
    """

}

workflow pre {
    take:
        pair_input_ch
    
    main:
        pair_input_ch.multiMap { sc_file, sp_file ->
                                sc_input: sc_file
                                sp_input: sp_file
                                }.set{ input }
        input.sc_input.view()
        input.sp_input.view()

}

/*
workflow {
  methods = "music,rctd,stereoscope"
  def python_methods = ["stereoscope"]
  // def methods_list = ["music", "rctd", "stereoscope"]
  methods_list = Arrays.asList(methods.split(','))
  //print(methods_list)
  print(!python_methods.disjoint(methods_list))
  
  // print(methods.disjoint(python_methods))
}
*/

workflow combine_ch {
  sc_input_ch = Channel.fromPath(params.sc_input)
  sp_input_ch = Channel.fromPath(params.sp_input)

  pair_input_ch = sc_input_ch.combine(sp_input_ch)
  pair_input_ch.view()
  //sc_input_ch.view()
  //sp_input_ch.view()
  // pre(pair_input_ch)
}

// A one-liner that I spent a very long time writing
bronze_standards = ((1..7)*8).sort().withIndex().collect{ it, i -> "bronze_standard_$it-${ -> i%8+1}".toString() }

all_pipelines = ["gold_standard_1", "gold_standard_2"] + bronze_standards

process bar {
  // container 'rocker/r-ver:3.6.3'
  echo true
  script:
    """
    # echo "task attempt: $task.attempt"
    # echo "cpus: $task.cpus"
    # echo \$(ls /mnt/d/)

    """
}
workflow testconfig {
  println(params.outdir.synvis)
}

workflow bleh{
  // bar()
  println(params.pipeline)
  // println(all_pipelines)
  // println(all_pipelines[3].getClass())
  if (!( params.pipeline ==~ /none/ )){
    println("hello world")
    test = "$params.rootdir/spotless-benchmark/data/reference/reference_${params.pipeline.split('-')[0]}.rds"
    //sc_input_ch = Channel.fromPath("$params.rootdir/spotless-benchmark/data/reference/reference_${params.pipeline.split('-')[0]}.rds")
    //sc_input_ch.view()
    println(test)
    if (!all_pipelines.contains(params.pipeline)){
      println("Pipeline not found.")

    }
    // sc_input = "$params.rootdir/spotless-benchmark/data/$params.pipeline/*[0-9].rds"
    // sp_input = "$params.rootdir/spotless-benchmark/data/$params.pipeline/reference*.rds"
    //sc_input_ch = Channel.fromPath(sc_input)
    //sc_input_ch.view()
    // sp_input_ch = Channel.fromPath(sp_input)
    // sp_input_ch.view()

    // sp_input_ch = Channel.fromPath(params.sp_input)
  }
}

synvis_types_map = [aud: "artificial_uniform_distinct", add: "artificial_diverse_distinct",
                    auo: "artificial_uniform_overlap", ado: "artificial_diverse_overlap",
                    adcd: "artificial_dominant_celltype_diverse", apdcd: "artificial_partially_dominant_celltype_diverse",
                    adrcd: "artificial_dominant_rare_celltype_diverse", arrcd: "artificial_regional_rare_celltype_diverse"]
synvis_types_fullnames = synvis_types_map.collect{ it.value }
synvis_types_flat = synvis_types_map.collect{[it.key, it.value]}.flatten()
/*
process foo2 {

}
*/
params.synvis = [type: "adrcd,blah,artificial_uniform_distinct,aud", reps: 10, n_regions: 5,
                region_var: "none", dataset_id: "1", n_spots_min: 50, n_spots_max: 500,
                visium_mean: 20000, visium_sd: 7000]

workflow {

  println("Hello world")
  // println(synvis_types)
  /*
  synvis_type_input = params.synvis.type.split(',').findAll{ synvis_types_flat.contains(it) }
                                        .collect{ ( synvis_types_fullnames.contains(it) ? it : synvis_types_map[it]) }
                                        .unique()
  println(params.synvis.findAll{ it.key != "type" && it.key != "reps" }.collect{ "--$it.key $it.value" }.join(" "))
  */
  // println(synvis_type_input)
   
  //println(synvis_type_input.collect{synvis_types[it]})
  // methods_list = Arrays.asList(methods.split(','))
  //println(synvis_types.collect{entry -> entry.key})
  //println(params.synvis_type)
}