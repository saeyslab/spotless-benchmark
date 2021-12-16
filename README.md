# SPAtial DEconvolution Benchmark
This is a repository for benchmarking spatial deconvolution and mapping tools.

You can either run the pipeline with the gold/silver/bronze standards or run it with your own data. In either case, you need to [install NextFlow](https://www.nextflow.io/docs/latest/getstarted.html), and either [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html).

In the future, you will be able to run this pipeline directly from NextFlow. From now you must clone the repository to run the pipeline.

‼**We highly recommend you to create a profile in the nextflow.config file**‼

There are currently three main profiles: 
- *local* uses the local environment 
- *prism* submits the job to a Sun Grid Engine cluster
- *hpc* submits the job to a high-performance computing cluster

When running locally, we suggest using it in conjuction with the *docker* profile, that is `nextflow run main.nf -profile local,docker`. You would need to modify `params.rootdir` (directory up to and including the repo, e.g., `"/home/$USER/spade-benchmark"`) and `workDir` (directory in which temporary files will be saved, you can also remove this) to suit your local folder structure.

If you want to deploy the pipeline in other clusters (e.g., an AWS cluster), you will have to create a new profile.

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
1. `generate_and_run` takes two single-cell Seurat objects, one to generate the synthetic data (`synvis.sc_input`) and one to use as input in deconvolution methods (`sc_input`). See the next section for more details.
2. `run_dataset` mode takes a single-cell Seurat object (`sc_input`), a directory containing the spatial dataset(s) (`sp_input`), and the cell type annotation column (`annot`). By default the spatial data is assumed to be generated using *synthvisium* and the annotation column is *celltype*. We can run the standards in this way also.

```
nextflow run main.nf -profile <profile_name> --mode run_dataset --sc_input standards/gold_standard_1/*.rds --sp_input standards/reference/gold_standard_1.rds --sp_type seqFISH
nextflow run main.nf -profile <profile_name> --mode run_dataset --sc_input standards/bronze_standard_1-1/*.rds --sp_input standards/reference/bronze_standard_1.rds
```

### Generating synthvisium data
The workflow `subworkflows/data_generation/generate_data.nf` generates synthetic visium data with *synthvisium*. The arguments are assumed to be stored in a dictionary, so it may be easier to provide this in a separate yaml/JSON file, shown below:

```
# synthvisium_params.yaml
synvis:
  sc_input: standards/reference/bronze_standard_1.rds
  clust_var: celltype
  reps: 3
  type: artificial_diverse_distinct,artificial_uniform_distinct
```
These parameters will return 6 synthetic datasets, with 3 replicates for each type. You can generate the data only (synthetic datasets will be copied to `outdir.synvis`), or run the whole pipeline immediately after.
```
# Only generate data
nextflow run subworkflows/data_generation/generate_data.nf -profile <profile_name> --params-file synthvisium_params.yaml

# Generate and run the whole pipeline
nextflow run main.nf -profile <profile_name> --mode generate_and_run --sc_input standards/reference/bronze_standard_1.rds --params-file synthvisium_params.yaml
```
In the second case, the same file will be used to generate synthetic data and to integrate with deconvolution methods. In our benchmark we use different files for this (akin to the training and test datasets in Machine Learning).

TODO: explain arguments?


## Platforms
The workflow has been tested on two platforms:
- NextFlow 21.04.3 on the Windows Subsystem for Linux (WSL2, Ubuntu 20.04), with Docker Desktop for Windows 4.1.0
- NextFlow 20.10.3 on CentOS Linux, with Singularity 3.8.1
