nextflow.enable.dsl=2

process convertRDStoH5AD {
    tag "convert_${rds_file_basename}"
    container 'csangara/seuratdisk:latest'
    label 'retry'

    input:
        path (rds_file)
    
    output:
        // Needs to return rds file for computing metrics
        tuple path ("${rds_file_basename}.h5ad"), path (rds_file)

    script:
        rds_file_basename = file(rds_file).getSimpleName()
        """
        Rscript $params.rootdir/subworkflows/deconvolution/convertRDStoH5AD.R \
        --input_path $rds_file
        """
}


process formatTSVFile {
    tag "format_${method_name}"
    container 'rocker/tidyverse:latest'
    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"

    input:
        tuple val(method_name), path (tsv_file), path(sp_input)
    output:
        tuple val(method_name), path (new_tsv_file), path(sp_input)
    script:
        new_tsv_file = file(tsv_file).getSimpleName()
        output_suffix = file(sp_input).getSimpleName()
        """
        #!/usr/bin/env Rscript
        deconv_matrix <- read.table("$tsv_file", sep="\t", header=TRUE, row.names=1)
        colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")
        deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]
        write.table(deconv_matrix, file="$new_tsv_file", sep="\t", quote=FALSE, row.names=FALSE)
        """

}