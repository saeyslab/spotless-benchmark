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

workflow {
  sc_input_ch = Channel.fromPath(params.sc_input)
  sp_input_ch = Channel.fromPath(params.sp_input)

  pair_input_ch = sc_input_ch.combine(sp_input_ch)
  pair_input_ch.view()
  //sc_input_ch.view()
  //sp_input_ch.view()
  // pre(pair_input_ch)
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
