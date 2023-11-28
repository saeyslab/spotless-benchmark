nextflow.enable.dsl=2
include {  fitSTRIDEModel } from './run_method.nf'
include { createDummyFile; formatTSVFile } from '../../helper_processes'

// This file is for running STRIDE on a pretrained model
// User must provide the zipped model file (--model_input) and txt file (--txt_input) from buildSTRIDEModel, in addition to the spatial data (--sp_input)
// Example use:
// nextflow run subworkflows/deconvolution/stride/fit_stride.nf -profile local,docker -params-file conf/method_params_silver.yaml \
// --model_input model.zip --txt_input txt_files.zip --sp_input "standards/silver_standard_1-1/*.h5ad"

workflow {

    // Print inputs (the timing isn't right with with view(), so do this instead)
    // Although view() has the advantage that it gives the absolute path (sc_input_ch.view())

    println("Model file:")
    println(file(params.model_input))

    println("\nTxt file:")
    println(file(params.txt_input))

    println("\nSpatial dataset(s):") 
    // With glob pattern, there will be multiple files
    params.sp_input =~ /\*/ ? file(params.sp_input).each{println "$it"} : println (file(params.sp_input))

    println("\nParameters:")
    println(params.deconv_args.stride)

    model_zip = Channel.fromPath(params.model_input)
    txt_zip = Channel.fromPath(params.txt_input)
    sp_input_ch = Channel.fromPath(params.sp_input) // Can have one or more files

    // Make model and txt_zip into a tuple
    model_txt_ch = model_zip.combine(txt_zip)

    // Double sp_input_ch to make a tuple
    sp_input_pair = createDummyFile(sp_input_ch)

    // Repeat the model for each sp_input
    model_txt_ch.combine(sp_input_pair)
        .multiMap { model_zip, txt_zip, sp_input, sp_input_rds ->
                        model: tuple model_zip, txt_zip
                        sp_input: tuple sp_input, sp_input_rds }
        .set{ stride_combined_ch }

    // Run STRIDE
    fitSTRIDEModel(stride_combined_ch.sp_input,
                   stride_combined_ch.model)

    formatTSVFile(fitSTRIDEModel.out) 

}
