process runSTRIDE {
    tag "stride_${sp_file_basename}"
    label "retry"
    label "longer_time"
    container 'csangara/sp_stride:latest'
    echo true

    input:
        tuple path (sc_input), path (sp_input), path (sp_input_rds)
    output:
        tuple val('stride'), path("$output"), path (sp_input_rds)

    script:
        sp_file_basename = file(sp_input).getSimpleName()
        output = "proportions_stride_${sp_file_basename}${params.runID_props}.preformat"
        args = ( params.deconv_args.stride ? params.deconv_args.stride : "" )

        """
        source activate stride
        echo "Creating annotation file..."
        python $params.rootdir/subworkflows/deconvolution/stride/createFiles.py \
           $sc_input $params.annot

        echo "Arguments: $args"
        if [[ "$args" == *"--markers"* ]] && [[ "$args" == *"--gene-use"* ]]
        then
          echo "ERROR: Both gene list and markers option is provided. You can only choose one."
          exit 1

        elif [[ "$args" == *"--markers"* ]]
        then
          echo "Identifying marker genes..."
          python $params.rootdir/subworkflows/deconvolution/stride/identify_markers.py \
            --sc-count $sc_input --annot $params.annot $args
          
          # Remove --markers from argument, add --gene-use
          re='(.*)--markers[ ]*[0-9]+(.*)'
          [[ "$args" =~ \$re ]]
          new_args="\${BASH_REMATCH[1]}\${BASH_REMATCH[2]} --gene-use markers.txt"

        else
          new_args="$args"
        fi

        STRIDE deconvolve --sc-count $sc_input --sc-celltype annot.txt \
            --st-count $sp_input \$new_args

        mv STRIDE_spot_celltype_frac.txt $output
        
        """
}

process buildSTRIDEModel {
    tag "stride_build_$tag_suffix"
    label "retry"
    label "longer_time"
    container 'csangara/sp_stride:latest'
    echo true

    input:
        path (sc_input)
    output:
        tuple path ("model.zip"), path("txt_files.zip")

    script:
        args = ( params.deconv_args.stride ? params.deconv_args.stride : "" )
        tag_suffix = file(sc_input).getSimpleName()
        
        """
        source activate stride
        echo "Creating annotation and dummy spatial files..."
        python $params.rootdir/subworkflows/deconvolution/stride/createFiles.py \
           $sc_input $params.annot dummy_arg

        STRIDE deconvolve --sc-count $sc_input --sc-celltype annot.txt \
            --st-count dummy_st.h5ad $args

        # The model/ folder contains multiple files, so zip it
        zip -r model.zip model
        zip txt_files.zip Model_selection.txt Gene_dict.txt
        """
}

process fitSTRIDEModel {
    tag "stride_$sp_file_basename"
    label "retry"
    container 'csangara/sp_stride:latest'
    echo true

    input:
        tuple path (sp_input), path (sp_input_rds)
        tuple path (model_zip), path(txt_files_zip)
    output:
        tuple val('stride'), path("$output"), path (sp_input_rds)
    script:
        sp_file_basename = file(sp_input).getSimpleName()
        output = "proportions_stride_${sp_file_basename}${params.runID_props}.preformat"
        args = ( params.deconv_args.stride ? params.deconv_args.stride : "" )
        
        """
        source activate stride
        unzip $model_zip
        unzip $txt_files_zip
        
        STRIDE deconvolve --model-dir model/ --st-count $sp_input $args
        
        mv STRIDE_spot_celltype_frac.txt $output
        """
}


