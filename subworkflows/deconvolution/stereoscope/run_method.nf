process buildStereoscopeModel {
    tag 'stereo_build'
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    container 'csangara/sp_stereoscope:latest'
    echo true

    input:
        // rds input actually is not needed
        tuple path (sc_input), path (sc_input_rds)
    output:
        tuple path ("R*.tsv"), path ("logits*.tsv")

    script:
        epochs = ( params.epoch_build ==~ /default/ ? "" : "-sce $params.epoch_build")
        args = ( params.deconv_args.stereoscope ? params.deconv_args.stereoscope : "" )
        gpu_flag = ( params.gpu ? "--gpu" : "" )
        println ("Building stereoscope model with ${ (params.gpu) ? "GPU" : "CPU" }...")

        """
        source activate stereoscope
        export CUDA_VISIBLE_DEVICES=$params.cuda_device

        stereoscope run --sc_cnt $sc_input --label_colname $params.annot \
            $epochs $args $gpu_flag -o \$PWD
        """
}

process fitStereoscopeModel {
    tag "stereo_$sp_file_basename"
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    container 'csangara/sp_stereoscope:latest'
    echo true

    input:
        tuple path (sp_input), path (sp_input_rds)
        tuple path (r_file), path (logits_file)
    output:
        tuple val('stereoscope'), path("$output"), path (sp_input_rds)
    script:
        sp_file_basename = file(sp_input).getSimpleName()
        output = "proportions_stereoscope_${sp_file_basename}${params.runID_props}.preformat"
        epochs = ( params.epoch_fit ==~ /default/ ? "" : "-ste $params.epoch_fit")
        args = ( params.deconv_args.stereoscope ? params.deconv_args.stereoscope : "" )
        gpu_flag = ( params.gpu ? "--gpu" : "" )

        println ("Received model files $r_file and $logits_file")
        println ("Fitting stereoscope model with ${ (params.gpu) ? "GPU" : "CPU" }...")
        
        """
        source activate stereoscope
        export CUDA_VISIBLE_DEVICES=$params.cuda_device
        
        stereoscope run --sc_fit $r_file $logits_file \
            --st_cnt $sp_input $epochs $args $gpu_flag -o \$PWD
        mv $sp_file_basename/W*.tsv $output
        """
}
