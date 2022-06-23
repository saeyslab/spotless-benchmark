process runDSTG {
    tag "dstg_$output_suffix"
    label "retry"
    container 'csangara/sp_dstg:latest'
    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"
    input:
        tuple path (sc_input), path (sp_input)

    output:
        tuple val('dstg'), path("$output"), path (sp_input)

    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_dstg_${output_suffix}${params.runID_props}"
        sc_filename = file(sc_input).getSimpleName()
        
        """
        WORKDIR=\$PWD
        
        # Convert Seurat object to matrix
        Rscript $params.rootdir/subworkflows/deconvolution/dstg/convertSeuratToMatrix.R \
            --sc_input $sc_input --sp_input $sp_input --annot $params.annot
        
        source activate dstg
        # Script doesn't run unless we're inside the directory
        cp -R ${params.deconv_args.dstg.dir}/DSTG/DSTG .
        cd DSTG/
        rm synthetic_data/*

        # Input: sc file, spatial file, and labels file
        Rscript convert_data.R ../${sc_filename}_matrix.rds \
        ../${output_suffix}_matrix.rds ../${sc_filename}_label.rds

        # Train the model
        python train.py

        # Add line numbers to output
        awk -v ln=1 '{print ln++  ","  \$0 }' DSTG_Result/predict_output.csv > tmp && mv tmp DSTG_Result/predict_output.csv
        
        # Add cell type labels
        { head -n 1 Datadir/Real_Label2.csv ; cat DSTG_Result/predict_output.csv; } > tmp && mv tmp DSTG_Result/predict_output.csv
        
        # Convert csv to tsv and move file to work directory
        sed -i -E 's/("([^"]*)")?,/\\2\t/g' DSTG_Result/predict_output.csv

        mv DSTG_Result/predict_output.csv ../$output
    
        """

}
