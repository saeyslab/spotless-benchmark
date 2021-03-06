# Spotless: A benchmark pipeline for <br> spatial deconvolution tools
This is a repository for running spatial deconvolution tools through a Nextflow pipeline. Currently, cell2location, DestVI, DSTG, MuSiC, NNLS, RCTD, SpatialDWLS, SPOTlight, stereoscope, STRIDE, Seurat, and Tangram have been implemented. To add your own method, see [Guideline_for_adding_deconvolution_tools.md](Guideline_for_adding_deconvolution_tools.md).

**Quickstart guide**
1. [Install NextFlow](https://www.nextflow.io/docs/latest/getstarted.html) and either [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html).
2. Clone this repository. (If Git LFS isn't installed, you will have to download `unit-test/test_sc_data.rds` and `unit-test/test_sp_data.rds` manually.)
3. Modify or create a profile in *nextflow.config*. To run the pipeline locally, modify `params.rootdir` under `profiles { local { ... } }` as the directory up to and including the reposity, e.g., `"/home/$USER/spotless-benchmark"`.
(The other two profiles are used inside computing clusters, namely *prism* for a Sun Grid Engine cluster, and *hpc* for a Slurm cluster.) 
4. While in the `spotless-benchmark/` directory:
```
nextflow run main.nf -profile local,docker --methods music --sc_input unit-test/test_sc_data.rds \
--sp_input unit-test/test_sp_data.rds --annot subclass

# If singularity is installed, use -profile local,singularity
```

This runs MuSiC as a test. The first run might take a few minutes because the containers have to be downloaded. If this works, you should see the proportions and metrics inside `deconv_proportions/` and `results/` respectively (these directories can be changed under `params.outdir`). 

To run more methods, type the method names separated with a comma but no spaces, e.g., `--methods rctd,music`. To adjust method parameters, see `subworkflows/deconvolution/README.md`.

## Running the pipeline
**Input:**
- Single-cell reference dataset: a Seurat (.rds) or Anndata (.h5ad) object containing cell type annotations in the object metadata
- Spatial dataset: a Seurat (.rds) object, Anndata (.h5ad) object, or named list of counts (see [Synthvisium object structure](#synthvisium-object-structure)) 

**Output:**
- A spot $\times$ cell type proportion matrix (tab-separated file) 
- Evaluation metrics (only if synthetic data follows the synthvisium object structure)

You can run the pipeline (`main.nf`) in three modes, `run_dataset` (the default mode, runs deconvolution tools on your own dataset), `generate_and_run` (generates synthetic datasets from your scRNA-seq data and benchmarks it), and `run_standard` (for reproducing our analysis with gold/bronze standards).
### *run_dataset*: running/benchmarking deconvolution tools on your own dataset
The `run_dataset` mode requires a single-cell object (`params.sc_input`), the path to the spatial dataset(s) (`params.sp_input`), and the cell type annotation column (`params.annot`, default: celltype). For real data, use the `skip_metrics` flag to skip the evaluation step.
```
nextflow run main.nf -profile <profile_name> --mode run_dataset --sc_input <PATH_TO_SC_FILE> --sp_input <PATH_TO_SP_FILE> \
--annot <ANNOT_COL> --methods <METHODS> --skip_metrics
```
We can also run the standards in this way.
```
nextflow run main.nf -profile <profile_name> --mode run_dataset --sp_input "standards/gold_standard_1/*.rds" \
--sc_input standards/reference/gold_standard_1.rds --annot celltype --methods <METHODS>
```
??? Don't forget to put any directories with [glob patterns](https://www.malikbrowne.com/blog/a-beginners-guide-glob-patterns) in quotes.

??? Metric calculation is only possible with an rds file following the [synthvisium object structure](#synthvisium-object-structure) (see `unit-test/test_sp_data.rds`).

### *generate_and_run*: generating and benchmarking your own synthetic datasets
`generate_and_run` takes two single-cell objects, one to generate the synthetic data (`params.synvis.sc_input`) and one to use as input in deconvolution methods (`params.sc_input`). It uses our *synthvisium* tool to generate synthetic data by running `subworkflows/data_generation/generate_data.nf`. Synthetic datasets will be copied to `params.outdir.synvis`. You can generate the data only or run the whole pipeline immediately after.
```
# Only generate data
nextflow run subworkflows/data_generation/generate_data.nf -profile <profile_name> -params-file conf/synthvisium.yaml

# Generate and run the whole pipeline
nextflow run main.nf -profile <profile_name> --mode generate_and_run --sc_input standards/reference/silver_standard_1.rds \
-params-file conf/synthvisium.yaml --methods <METHODS>
```
The arguments to synthvisium are best provided in a separate yaml/JSON file. Check out `conf/synthvisium.yaml` for a detailed description of arguments. Minimally, you need four arguments:
```
# conf/synthvisium.yaml
synvis:
  sc_input: standards/reference/silver_standard_1.rds             # single-cell reference input
  clust_var: celltype                                             # name of metadata column with cell type annotation
  reps: 3                                                         # number of replicates per dataset type (abundance pattern)
  type: artificial_diverse_distinct,artificial_uniform_distinct   # dataset types
```
This will return 3 replicates for each dataset type, resulting in 6 files. You can also adjust other parameters such as the number of spots and mean or standard deviation per spot. Note that in this example, the same file was used to generate synthetic data and to integrate with deconvolution methods. In our benchmark we use different files for this (akin to the training and test datasets in machine learning).

### *run_standard*: reproducing our analysis
Download the datasets from Zenodo and extract the file.
```
cd spotless-benchmark
wget https://zenodo.org/record/6786528/files/standards.tar.gz?download=1 -O standards.tar.gz
tar -xvzf standards.tar.gz 
```
Then run the pipeline with the `run_standard` mode.
```
nextflow run main.nf -profile <profile_name> --mode run_standard --standard <standard_name> -c standards/standard.config \
--methods <METHODS>
```
All folder names (except `reference`) can be used as the *standard_name*. For instance, to run the gold standard of seqFISH+ cortex dataset or the brain cortex silver standard, you would do
```
nextflow run main.nf -profile <profile_name> --mode run_standard --standard gold_standard_1 -c standards/standard.config
nextflow run main.nf -profile <profile_name> --mode run_standard --standard silver_standard_1-1 -c standards/standard.config
```

## Pipeline arguments (Advanced use)
You can find the default arguments of the pipeline in the `nextflow.config` file, under the `params` scope. These can be overwritten by parameters provided in the command line or in an external JSON/YAML file (see exact priorities [here](https://www.nextflow.io/docs/latest/config.html)).
* `methods`: deconvolution methods to run in the pipeline, must be comma-separated with no space, e.g.,  <br /> `--methods music,rctd` (case-insensitive, default: all)
* `mode`: as explained above, the different modes in which the pipeline can be run (run_standard, run_dataset, generate_and_run)
* `annot`: the cell type annotation column in the input scRNA-seq Seurat object (default: celltype)
* `outdir`: location to save the proportions, metrics, and synthetic data (default: `deconv_proportions/`, `results/`, `synthetic_data/`). Best to define this under your profiles.
* `sampleID`: the column containing batch information for the input scRNA-seq Seurat object (default: none) 
* `deconv_args`: extra parameters to pass onto deconvolution algorithms (default: []). For a more detailed explanation and list of parameters for each method, see `subworkflows/deconvolution/README.md`.
* `synvis`: synthvisium arguments, see `subworkflows/synthvisium.yaml`
* `gpu`: add this flag to use host GPU, see below
* `verbose`: add this flag to print input files
* `skip_metrics`: add this flag to skip the final step which calculates evaluation metrics
* `runID_props`, `runID_metrics`: a suffix added to the proportions and metrics file output by the pipeline
* `remap_annot`: if you want to use another celltype annotation to calculate the metrics, see below

### GPU usage
Stereoscope, cell2location, DestVI, and Tangram can make use of a GPU to shorten their runtimes. You can do this by providing the `--gpu` flag when running the pipeline. If you have a specific GPU you want to use, you will have to provide the index with `cuda_device` (default: 0). This works from inside the containers, so you still do not have to install the programs locally. However, the containers were built on top of a [NVIDIA base image](https://hub.docker.com/r/nvidia/cuda) with CUDA version 10.2, so your GPU must be compatible with this CUDA version.
```
nextflow run main.nf -profile <profile_name> --mode run_standard --standard gold_standard_1 \
-c standards/standard.config --methods stereoscope,cell2location --gpu #--cuda_device 1
```

## Running selected parts of the pipeline
It is also possible to only run a subworkflow instead of the whole pipeline.

### Generating synthvisium data
See __*generate_and_run*: generating and benchmarking your own synthetic datasets__.

Briefly, running the following code will save an rds file to `params.outdir.synvis`:
```
nextflow run subworkflows/data_generation/generate_data.nf -profile <profile_name> -params-file conf/synthvisium.yaml
```

#### Synthvisium object structure
The output of synthvisium is a named list of matrices and lists. There are three necessary components:
1) **counts**: a gene $\times$ spot count matrix (preferably raw counts)
2) **relative_spot_composition**: a spot $\times$ cell types relative proportion matrix

You can look at any of the files in `standards/silver_standard` for a better idea.

### Running methods
```
nextflow run subworkflows/deconvolution/run_methods.nf -profile <profile_name> \
--sc_input <scRNAseq_dataset> --sp_input <spatial_dataset> --annot <celltype annotation column> \
--methods <METHODS>
```
By default, all methods are run (equivalent to giving an argument `--methods music,rctd,spatialdwls,spotlight,stereoscope,cell2location,destvi,dstg,tangram,seurat,stride`). To run only certain methods, make sure the values are comma-separated without any spaces.

### Computing metrics
It is possible to add more metrics in the `subworkflows/evaluation/metrics.R` yourself, then recalculate the metrics without running the entire pipeline. You have to provide the ground truth datasets with the synthvisium object structure (`params.sp_input`), and the script will look for the proportions at `params.outdir.props`.
```
nextflow run subworkflows/evaluation/evaluate_methods.nf -profile <profile_name> \
  --sp_input "standards/gold_standard_1/*.rds"
```

### Remapping cell type annotations
You can sum up the proportions of closely related cell types together to see if the methods perform better with a broader annotation. You can do this by creating a tsv file that maps one cell type to another (see `standards/gold_standard_1/conversion.tsv` for an example). The deconvolution will still be performed with the original annotations, but the proportions will be summed during metrics calculation. You have to provide the absolute file path with the `remap_annot` parameter. 
```
nextflow run main.nf -profile <profile_name> --mode run_standard --standard gold_standard_1 \
-c standards/standard.config --remap_annot /home/$USER/spotless-benchmark/standards/gold_standard_1/conversion.tsv

# Only rerun the calculations - add file suffix to not overwrite existing metrics file
nextflow run subworkflows/evaluation/evaluate_methods.nf -profile <profile_name> \
  --sp_input "standards/gold_standard_1/*.rds" --runID_metrics "_coarse" \
  --remap_annot /home/$USER/spotless-benchmark/standards/gold_standard_1/conversion.tsv 
```

## Platforms
The workflow has been tested on three platforms:
- NextFlow 21.04.3 on the Windows Subsystem for Linux (WSL2, Ubuntu 20.04), with Docker Desktop for Windows 4.1.0
- NextFlow 20.10.3 on CentOS 7, with Singularity 3.8.1
- NextFlow 21.03.0 on CentOS 7, with Singularity 3.8.5 (and NVIDIA Volta V100 GPUs)
