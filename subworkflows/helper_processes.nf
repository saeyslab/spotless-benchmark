nextflow.enable.dsl=2

process convertBetweenRDSandH5AD {
    tag "convert_${file_basename}"
    container 'csangara/seuratdisk:latest'
    label 'retry'

    input:
        path (file_to_convert)
    
    output:
        // Needs to return rds file for computing metrics
        tuple path ("${file_basename}.h5ad"), path ("${file_basename}.rds")

    script:
        file_basename = file(file_to_convert).getSimpleName()
        """
        Rscript $params.rootdir/subworkflows/deconvolution/convertBetweenRDSandH5AD.R \
        --input_path $file_to_convert
        """
}


process formatTSVFile {
    tag "format_${method_name}"
    label 'trivial'
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

process createDummyFile {
    tag "dummy_${file_name}"
    label 'trivial'

    input:
        path (h5ad_file)
    output:
        tuple path (h5ad_file), path ("${file_name}.rds")

    script:
        file_name = file(h5ad_file).getSimpleName()
        """
        touch ${file_name}.rds
        """
}

// This is for if you want to convert files outside the normal workflow
process convertProcessForWorkflow {
    tag "convert_${file_basename}"
    container 'csangara/seuratdisk:latest'
    publishDir { "${file_directory}/" },
               mode: 'copy', pattern: "*.${ext_to_include}"
    label 'retry'

    input:
        path (file_to_convert)
    
    output:
        // Needs to return rds file for computing metrics
        tuple path ("${file_basename}.h5ad"), path ("${file_basename}${output_file_ext}.rds")

    script:
        file_basename = file(file_to_convert).getSimpleName()
        file_directory = params.convert_input =~ /\*/ ? file(params.convert_input)[0].getParent() :
                                                        file(params.convert_input).getParent()
        ext_to_include = file(file_to_convert).getExtension().toLowerCase() =~ /h5/ ? "rds" : "h5ad"
        output_file_ext = ext_to_include ==~ /rds/ ? "_fromH5AD" : ""
        println("Saving ${ext_to_include} file to ${file_directory}")
        """
        Rscript $params.rootdir/subworkflows/deconvolution/convertBetweenRDSandH5AD.R \
        --input_path $file_to_convert --output_file_ext $output_file_ext
        """
}

workflow convertWorkflow {
    println("Input files:")
    params.convert_input =~ /\*/ ? file(params.convert_input).each{println "$it"} : println (file(params.convert_input))
    input_files = Channel.fromPath(params.convert_input)
    convertProcessForWorkflow(input_files)

}
