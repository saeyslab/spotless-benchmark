We implemented two ways to run STRIDE in `run_method.nf`:
1. Model building + model fitting (like cell2location or stereoscope) -> `buildSTRIDEModel`, `fitSTRIDEModel`
2. Both steps at once -> `runSTRIDE`

By default, we used the latter implementation because it allows for running STRIDE with only marker genes (implemented by us with the parameter `--markers`). When the model building step takes really long however, it may be favorable to use the former implementation instead. This must be done manually in `../run_methods.nf`, under the STRIDE code block:

```
 if ( methods =~ /stride/ ) {
                
                /*
                buildSTRIDEModel(sc_input_conv)

                // Repeat model output for each spatial file
                buildSTRIDEModel.out.combine(sp_input_pair)
                .multiMap { zip_file, txt_files, sp_file_h5ad, sp_file_rds ->
                            model: tuple zip_file, txt_files
                            sp_input: tuple sp_file_h5ad, sp_file_rds }
                .set{ stride_combined_ch }
                
                fitSTRIDEModel(stride_combined_ch.sp_input,
                               stride_combined_ch.model)
                formatSTRIDE(fitSTRIDEModel.out) 
                */
                
                
                runSTRIDE(sc_input_conv.combine(sp_input_pair))
                formatSTRIDE(runSTRIDE.out)
                
                output_ch = output_ch.mix(formatSTRIDE.out)
            }
```
Users can uncomment the commented section, and instead comment the two lines containing `runSTRIDE` and `formatSTRIDE`.

We also implemented an extra workflow in `fit_stride.nf`, which allows you to run STRIDE on a pretrained model. You must provide the path to the model.zip and txt_files.zip files returned from `buildSTRIDEModel` with the arguments `--model_input` and `--txt_input`, respectively. An example use case is as follows:
```
nextflow run subworkflows/deconvolution/stride/fit_stride.nf -profile local,docker -params-file conf/method_params_silver.yaml \
--model_input model.zip --txt_input txt_files.zip --sp_input "standards/silver_standard_1-1/*.h5ad"
```

This workflow only supports h5ad files as the spatial input and does not compute any evaluation metrics.
