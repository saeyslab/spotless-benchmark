# Spotless: A benchmark pipeline for <br> spatial deconvolution tools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8211492.svg)](https://doi.org/10.5281/zenodo.8211492)

This is a repository for running spatial deconvolution tools through a Nextflow pipeline, as described in [Sang-aram et al. (2024)](https://elifesciences.org/articles/88431). <br> Currently, cell2location, DestVI, DSTG, MuSiC, NNLS, RCTD, SpatialDWLS, SPOTlight, stereoscope, STRIDE, Seurat, and Tangram have been implemented. To add your own method, see [Guideline_for_adding_deconvolution_tools.md](Guideline_for_adding_deconvolution_tools.md).

**Quickstart guide**
1. [Install NextFlow](https://www.nextflow.io/docs/latest/getstarted.html) and either [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html).


2. Download [containers](https://hub.docker.com/u/csangara).

<details> 
<summary>Sample code</summary>


```
# Method containers
for method in cell2location destvi dstg music nnls rctd spatialdwls spotlight stereoscope stride seurat tangram
do
docker pull csangara/sp_${method}:latest
done

# Other
docker pull csangara/sp_eval:latest # for computing performance metrics
docker pull csangara/seuratdisk:latest # for converting between Seurat and AnnData objects internally

# For singularity
singularity pull docker://csangara/${container_name}:latest
# You might have to define path to singularity containers using the "singularity.cacheDir" directive in the config file
```


</details>

3. Clone this repository. Check locally that `unit-test/test_sc_data.rds` is 82.7 MB and `unit-test/test_sp_data.rds` is 290 KB. If not, manually download them [[1]](https://github.com/saeyslab/spotless-benchmark/raw/main/unit-test/test_sc_data.rds)[[2]](https://github.com/saeyslab/spotless-benchmark/raw/main/unit-test/test_sp_data.rds) into the `unit-test` folder. (This occurs when Git LFS is not installed.)
5. Modify or create a profile in *nextflow.config*. To run the pipeline locally, modify `params.rootdir` under `profiles { local { ... } }` as the directory up to and including the reposity, e.g., `"/home/$USER/spotless-benchmark"`.
(The other two profiles are used inside computing clusters: *prism* for a Sun Grid Engine cluster, and *hpc* for a Slurm cluster.) 
5. While in the `spotless-benchmark/` directory:
```
nextflow run main.nf -profile local,docker --methods music --sc_input unit-test/test_sc_data.rds \
--sp_input unit-test/test_sp_data.rds --annot subclass

# If singularity is installed, use -profile local,singularity
```

This runs MuSiC as a test. If this works, you should see the proportions and metrics inside `deconv_proportions/` and `results/` respectively (these directories can be changed under `params.outdir`).

To run more methods, type the method names separated with a comma but no spaces, e.g., `--methods rctd,music`. To adjust method parameters, see `subworkflows/deconvolution/README.md`.

‼**WARNING:** Do not run multiple methods simultaneously on your local computer, please only do so on a cluster!

## Running the pipeline
**Input:**
- Single-cell reference dataset: a Seurat (.rds) or Anndata (.h5ad) object containing cell type annotations in the object metadata
- Spatial dataset: a Seurat (.rds) object, Anndata (.h5ad) object, or named list of counts (see [Synthspot object structure](#synthspot-object-structure)) 

**Output:**
- A spot $\times$ cell type proportion matrix (tab-separated file) 
- Evaluation metrics (only if synthetic data follows the synthspot object structure)

You can run the pipeline (`main.nf`) in three modes, `run_dataset` (the default mode, runs deconvolution tools on your own dataset), `generate_and_run` (generates synthetic datasets from your scRNA-seq data and benchmarks it), and `run_standard` (for reproducing our analysis with gold/silver standards).
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
‼ Don't forget to put any directories with [glob patterns](https://www.malikbrowne.com/blog/a-beginners-guide-glob-patterns) in quotes.

‼ Metric calculation is only possible with an rds file following the [Synthspot object structure](#synthspot-object-structure) (see `unit-test/test_sp_data.rds`).

### *generate_and_run*: generating and benchmarking your own synthetic datasets
`generate_and_run` takes two single-cell objects, one to generate the synthetic data (`params.synthspot.sc_input`) and one to use as input in deconvolution methods (`params.sc_input`). It uses our *synthspot* tool to generate synthetic data by running `subworkflows/data_generation/generate_data.nf`. Synthetic datasets will be copied to `params.outdir.synthspot`. You can generate the data only or run the whole pipeline immediately after.
```
# Only generate data
nextflow run subworkflows/data_generation/generate_data.nf -profile <profile_name> -params-file conf/synthspot.yaml

# Generate and run the whole pipeline
nextflow run main.nf -profile <profile_name> --mode generate_and_run --sc_input standards/reference/silver_standard_1.rds \
-params-file conf/synthspot.yaml --methods <METHODS>
```
The arguments to synthspot are best provided in a separate yaml/JSON file. Check out `conf/synthspot.yaml` for a detailed description of arguments. Minimally, you need four arguments:
```
# conf/synthspot.yaml
synthspot:
  sc_input: standards/reference/silver_standard_1.rds             # single-cell reference input
  clust_var: celltype                                             # name of metadata column with cell type annotation
  reps: 3                                                         # number of replicates per dataset type (abundance pattern)
  type: artificial_diverse_distinct,artificial_uniform_distinct   # dataset types
```
This will return 3 replicates for each dataset type, resulting in 6 files. You can also adjust other parameters such as the number of spots and mean or standard deviation per spot. Note that in this example, the same file was used to generate synthetic data and to integrate with deconvolution methods. In our benchmark we use different files for this (akin to the training and test datasets in machine learning).

### *run_standard*: reproducing our analysis
Download and extract `standards.tar.gz` from [Zenodo](https://zenodo.org/record/8211492).
```
cd spotless-benchmark
wget https://zenodo.org/record/8211492/files/standards.tar.gz?download=1 -O standards.tar.gz
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

### Running the liver dataset
Download and extract `liver_dataset.tar.gz` from [Zenodo](https://zenodo.org/record/8211492). `liver_README.txt` contains a description of the files. You could run the analysis by using `--mode run_dataset`, but we also provide the `conf/liver_mouse_visium.config` file which contains the parameters we used for each method as well as an additional argument (`ref_type`) you could add to automatically select a certain reference dataset. (You will however have to replace the base directories of `sc_input` and `sp_input` in this config file.) The following code snippet uses the Nuc-seq dataset for deconvolution:

```
nextflow run main.nf -profile hpc --mode run_dataset --methods rctd -c conf/liver_mouse_visium.config \
--ref_type nuclei  # options: nuclei, inVivo, exVivo, noEC, 9celltypes \
--annot annot_cd45 # options: annot, annot_fine, annot_cd45 \
--file_type rds    # options: rds, h5ad \
#--gpu
```
Note that it is much more efficient to directly use the appropriate file types (rds/h5ad) rather than using the internal conversion. Therefore, it is best to run R and Python methods separately.

## Pipeline arguments (Advanced use)
You can find the default arguments of the pipeline in the `nextflow.config` file, under the `params` scope. These can be overwritten by parameters provided in the command line or in an external JSON/YAML file (see exact priorities [here](https://www.nextflow.io/docs/latest/config.html)).
* `methods`: deconvolution methods to run in the pipeline, must be comma-separated with no space, e.g.,  <br /> `--methods music,rctd` (case-insensitive, default: all)
* `mode`: as explained above, the different modes in which the pipeline can be run (run_standard, run_dataset, generate_and_run)
* `annot`: the cell type annotation column in the input scRNA-seq Seurat object (default: celltype)
* `outdir`: location to save the proportions, metrics, and synthetic data (default: `deconv_proportions/`, `results/`, `synthetic_data/`). Best to define this under your profiles.
* `sampleID`: the column containing batch information for the input scRNA-seq Seurat object (default: none) 
* `deconv_args`: extra parameters to pass onto deconvolution algorithms (default: []). For a more detailed explanation and list of parameters for each method, see `subworkflows/deconvolution/README.md`.
* `synthspot`: synthspot arguments, see `subworkflows/synthspot.yaml`
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

### Generating synthspot data
See __*generate_and_run*: generating and benchmarking your own synthetic datasets__.

Briefly, running the following code will save an rds file to `params.outdir.synthspot`:
```
nextflow run subworkflows/data_generation/generate_data.nf -profile <profile_name> -params-file conf/synthspot.yaml
```

#### Synthspot object structure
The output of synthspot is a named list of matrices and lists. There are three necessary components:
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
It is possible to add more metrics in the `subworkflows/evaluation/metrics.R` yourself, then recalculate the metrics without running the entire pipeline. You have to provide the ground truth datasets with the synthspot object structure (`params.sp_input`), and the script will look for the proportions at `params.outdir.props`.
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

### Converting between Seurat/AnnData (rds/h5ad) objects
The workflow detects whether the input file is of rds or h5ad format and converts the file to its counterpart. Use `params.convert_input` to indicate the  file(s) to be converted. In case of multiple files, they cannot be from different folders.
```
nextflow run subworkflows/helper_processes.nf -entry convertWorkflow -profile <profile_name> \
--convert_input "standards/silver_standard_2-1/*.rds"
```

## Platforms
The workflow has been tested on the following platforms:
- Local: NextFlow 21.04.3 on Windows Subsystem for Linux (WSL2, Ubuntu 20.04), with Docker Desktop for Windows 4.1.0
- Local: Nextflow 21.10.6 on CentOS 8, with Docker 20.10.14
- SGE cluster: NextFlow 20.10.3 on CentOS 7, with Singularity 3.8.1 
- Slurm cluster: NextFlow 21.03.0 on CentOS 7, with Singularity 3.8.5 (and NVIDIA Volta V100 GPUs)

## Citation
If you used our pipeline in your publication, please cite:

Sang-Aram, C., Browaeys, R., Seurinck, R., & Saeys, Y. (2024). Spotless, a reproducible pipeline for benchmarking cell type deconvolution in spatial transcriptomics. Elife, 12, RP88431. https://doi.org/10.7554/eLife.88431.3
