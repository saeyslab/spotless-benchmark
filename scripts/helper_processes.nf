nextflow.enable.dsl=2

process convertRDStoH5AD {
    container 'csangara/seuratdisk:latest'

    input:
        path (rds_file)
        val (input_type)
    
    output:
        file ("${rds_file_basename}.h5ad")

    script:
        rds_file_basename = file(rds_file).getSimpleName()
        """
        Rscript $params.rootdir/scripts/deconvolution/convertRDStoH5AD.R \
        --input_path $rds_file --input_type $input_type
        """
}


process formatTSVFile {
    publishDir params.outdir.props, mode: 'copy'

    input:
        tuple val(method_name), path (tsv_file)
    output:
        tuple val(method_name), path (new_tsv_file)
    script:
        new_tsv_file = file(tsv_file).getSimpleName()
        """
        #!/usr/bin/env Rscript
        deconv_matrix <- read.table("$tsv_file", sep="\t", header=TRUE, row.names=1)
        colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")
        deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]
        write.table(deconv_matrix, file="$new_tsv_file", sep="\t", quote=FALSE, row.names=FALSE)
        """

}