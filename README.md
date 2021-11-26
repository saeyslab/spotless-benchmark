# SPAtial DEconvolution Benchmark
This is a repository for benchmarking spatial deconvolution and mapping tools.

You can either run the pipeline with the gold/silver/bronze standards or run it with your own data. In either case, you need to [install NextFlow](https://www.nextflow.io/docs/latest/getstarted.html), and either [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html).

In the future, you will be able to run this pipeline directly from NextFlow. From now you must clone the repository to run the pipeline.

‼**We highly recommend you to create a profile in the nextflow.config file**‼

There are currently three profiles: 
- *local_env* uses the local environment without docker containers (i.e., you have to install all methods locally)
- *local_docker* uses the local environment with docker containers
- *prism* submits the job to a Sun Grid Engine cluster

When running locally, we suggest using *local_docker* by modifying `params.rootdir` (directory preceding the cloned repository) and `workDir` (directory in which temporary files will be saved, you can also remove this).

If you want to deploy the pipeline in other clusters (e.g., HPC or AWS cluster), you will have to create a new profile.

## Reproducing the pipeline with standards
First, download the datasets from Zenodo and place them in the `data/` folder.
```
cd spade-benchmark
wget <link>
tar -xf standards.tar.gz -C data/
```
To run the pipeline on a dataset:
```
nextflow run main.nf -profile <profile_name> --mode run_pipeline --pipeline <pipeline_name> -c data/pipeline.config
```
All folder names (except `reference`) can be used as the pipeline_name. For instance, to run the gold standard of seqFISH+ cortex dataset or the brain cortex bronze standard, you would do
```
nextflow run main.nf -profile <profile_name> --mode run_pipeline --pipeline gold_standard_1 -c data/pipeline.config
nextflow run main.nf -profile <profile_name> --mode run_pipeline --pipeline bronze_standard_1-1 -c data/pipeline.config
```

## Running the pipeline on your own dataset
Run the pipeline with the `run_dataset` mode. At the minimum you would need to provide a single-cell Seurat object (`sc_input`), a directory containing the spatial dataset(s) (`sp_input`), and the cell type annotation column (`annot`). By default the spatial data is assumed to be generated using *synthvisium* and the annotation column is *celltype*. We can run the standards in this way also.

```
nextflow run main.nf -profile <profile_name> --mode run_dataset --sc_input data/gold_standard_1/*.rds --sp_input data/reference/gold_standard_1.rds --sp_type seqFISH
nextflow run main.nf -profile <profile_name> --mode run_dataset --sc_input data/bronze_standard_1-1/*.rds --sp_input data/reference/bronze_standard_1.rds
```

### Generating synthvisium data
TODO


## Platforms
The workflow has been tested on two platforms:
- NextFlow 21.04.3 on the Windows Subsystem for Linux (WSL2, Ubuntu 20.04), with Docker Desktop for Windows 4.1.0
- NextFlow 20.10.3 on CentOS Linux, with Singularity 3.8.1
