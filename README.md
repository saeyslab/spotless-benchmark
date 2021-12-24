# SPAtial DEconvolution Benchmark
This is a repository for benchmarking spatial deconvolution and mapping tools.

You can either run the pipeline with the gold/silver/bronze standards or run it with your own data. In either case, you need to [install NextFlow](https://www.nextflow.io/docs/latest/getstarted.html), and either [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html).

In the future, you will be able to run this pipeline directly from NextFlow. For now you must clone the repository to run the pipeline.

‼**You need to modify/create a profile to run the pipeline**‼
Currently, there are three profiles: *local* uses the local environment to execute the pipeline, *prism* uses a Sun Grid Engine cluster, and *hpc* uses a Slurm cluster.

To run the pipeline locally, modify `params.rootdir` as the directory up to and including the reposity, e.g., `"/home/$USER/spade-benchmark"`. To use containers, you also need to modify the bind mounts of `docker.runOptions`. Then you can run the pipeline with `nextflow run main.nf -profile local,docker`.

## Reproducing the pipeline with standards
First, download the datasets from Zenodo and place them in the `standards/` folder. The file is 5GB.
```
cd spade-benchmark
wget https://zenodo.org/record/5727614/files/standards.tar.gz?download=1 -O standards.tar.gz
tar -xf standards.tar.gz -C standards/
mv standards/data/* standards/ ; rmdir standards/data # file structure to be fixed in the future
rm standards.tar.gz
```
Then run the pipeline with the `run_standard` mode.
```
nextflow run main.nf -profile <profile_name> --mode run_standard --standard <standard_name> -c standards/standard.config
```
All folder names (except `reference`) can be used as the *standard_name*. For instance, to run the gold standard of seqFISH+ cortex dataset or the brain cortex bronze standard, you would do
```
nextflow run main.nf -profile <profile_name> --mode run_standard --standard gold_standard_1 -c standards/standard.config
nextflow run main.nf -profile <profile_name> --mode run_standard --standard bronze_standard_1-1 -c standards/standard.config
```

## Running the pipeline on your own dataset
There are two modes:
1. `run_dataset` mode takes a single-cell Seurat object (`sc_input`), a directory containing the spatial dataset(s) (`sp_input`), and the cell type annotation column (`annot`). By default the spatial data is assumed to be generated using *synthvisium* and the annotation column is *celltype*. We can run the standards in this way also.
```
nextflow run main.nf -profile <profile_name> --mode run_dataset --sp_input "standards/gold_standard_1/*.rds" \
--sc_input standards/reference/gold_standard_1.rds --sp_type seqFISH

nextflow run main.nf -profile <profile_name> --mode run_dataset --sp_input "standards/bronze_standard_1-1/*.rds" \
--sc_input standards/reference/bronze_standard_1_brain_cortex.rds
```
‼ Don't forget to put any directories with [glob patterns](https://www.malikbrowne.com/blog/a-beginners-guide-glob-patterns) in quotes.

2. `generate_and_run` takes two single-cell Seurat objects, one to generate the synthetic data (`synvis.sc_input`) and one to use as input in deconvolution methods (`sc_input`). See the next section for more details.

### Generating synthvisium data
The workflow `subworkflows/data_generation/generate_data.nf` generates synthetic visium data with *synthvisium*. The arguments are assumed to be stored in a dictionary, so it may be easier to provide this in a separate yaml/JSON file. The four required arguments are:
```
# synthvisium_params.yaml
synvis:
  sc_input: standards/reference/bronze_standard_1.rds
  clust_var: celltype
  reps: 3
  type: artificial_diverse_distinct,artificial_uniform_distinct
```
These parameters will return 3 replicates for each dataset type, resulting in 6 files. You can also adjust other parameters such as the number of spots and mean or standard deviation per spot (see `subworkflows/data_generation/generate_synthetic_data.R`). You can generate the data only (synthetic datasets will be copied to `outdir.synvis`), or run the whole pipeline immediately after.
```
# Only generate data
nextflow run subworkflows/data_generation/generate_data.nf -profile <profile_name> --params-file synthvisium_params.yaml

# Generate and run the whole pipeline
nextflow run main.nf -profile <profile_name> --mode generate_and_run --sc_input standards/reference/bronze_standard_1.rds --params-file synthvisium_params.yaml
```
In the second case, the same file will be used to generate synthetic data and to integrate with deconvolution methods. In our benchmark we use different files for this (akin to the training and test datasets in machine learning).

## Pipeline arguments (Advanced use)
You can find the default arguments of the pipeline in the `nextflow.config` file, under the `params` scope. These can be overwritten by parameters provided in the command line or in an external JSON/YAML file (see exact priorities [here](https://www.nextflow.io/docs/latest/config.html)).
* `methods`: deconvolution methods to run in the pipeline, must be in small laters and comma-separated with no space, e.g.,  <br /> `--methods music,rctd` (default: all)
* `mode`: as explained above, the different modes in which the pipeline can be run (run_standard, run_dataset, generate_and_run)
* `annot`: the cell type annotation column in the input scRNA-seq Seurat object (default: subclass)
* `sampleID`: the column containing batch information for the input scRNA-seq Seurat object (default: none) 
* `deconv_args`: extra parameters to pass onto deconvolution algorithms (default: []). For a syntax example, check out `conf/test.config`. Can also be passed with the command line, e.g., `--deconv_args.cell2location "-p 10"`
* `synvis`: synthvisium arguments, see "Generating synthvisium data"
* `gpu`: add this flag to use host GPU, see below
* `verbose`: add this flag to print input files
* `runID_props`, `runID_metrics`: a suffix added to the proportions and metrics file output by the pipeline
* `remap_annot`: if you want to use another celltype annotation to calculate the metrics, see below

### GPU usage
Stereoscope and cell2location can make use of a GPU to shorten their runtimes. You can do this by providing the `--gpu` flag when running the pipeline. If you have a specific GPU you want to use, you will have to provide the index with `cuda_device` (default: 0). This works from inside the containers, so you still do not have to install the programs locally (the containers were built on top of a [NVIDIA base image](https://hub.docker.com/r/nvidia/cuda)).
```
nextflow run main.nf -profile <profile_name> --mode run_standard --standard gold_standard_1 \
-c standards/standard.config --methods stereoscope,cell2location --gpu #--cuda_device 1
```

### Remapping cell type annotations
You can sum up the proportions of closely related cell types together to see if the methods perform better with a broader annotation. You can do this by creating a tsv file that maps one cell type to another (see `standards/gold_standard_1/conversion.tsv` for an example). The deconvolution will still be performed with the original annotations, but the proportions will be summed during metrics calculation. You have to provide the absolute file path with the `remap_annot` parameter. 
```
nextflow run main.nf -profile <profile_name> --mode run_standard --standard gold_standard_1 \
-c standards/standard.config --remap_annot /home/$USER/spade-benchmark/standards/gold_standard_1/conversion.tsv

# Only rerun the calculations - add file suffix to not overwrite existing metrics file
nextflow run subworkflows/evaluation/evaluate_methods.nf -profile <profile_name> \
  --sp_input "standards/gold_standard_1/*.rds" --sp_type seqFISH --runID_metrics "_coarse" \
  --remap_annot /home/$USER/spade-benchmark/standards/gold_standard_1/conversion.tsv 
```

## Platforms
The workflow has been tested on three platforms:
- NextFlow 21.04.3 on the Windows Subsystem for Linux (WSL2, Ubuntu 20.04), with Docker Desktop for Windows 4.1.0
- NextFlow 20.10.3 on CentOS 7, with Singularity 3.8.1
- NextFlow 21.03.0 on CentOS 7, with Singularity 3.8.5 (and NVIDIA Volta V100 GPUs)
