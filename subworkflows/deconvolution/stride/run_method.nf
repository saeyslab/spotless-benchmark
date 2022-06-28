process buildSTRIDEModel {
    tag 'stride_build'
    label "retry"
    container 'csangara/sp_stride:latest'
    echo true

    input:
        path (sc_input)
    output:
        tuple path ("model.zip"), path("txt_files.zip")

    script:
        args = ( params.deconv_args.stride ? params.deconv_args.stride : "" )

        """
        source activate stride
        echo "Creating annotation and dummy spatial files..."
        python $params.rootdir/subworkflows/deconvolution/stride/createFiles.py \
           $sc_input $params.annot

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
