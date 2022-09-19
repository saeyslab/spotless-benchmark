## Method parameters

- [cell2location](#cell2location)
- [DestVI](#destvi)
- [MuSiC](#music)
- [Seurat](#seurat)
- [spatialDWLS](#spatialdwls)
- [SPOTlight](#spotlight)
- [stereoscope](#stereoscope)
- [STRIDE](#stride)
- [Tangram](#tangram)

(RCTD and DSTG have no adjustable parameters.)

### Usage 
Parameters for each method can be provided in a *yaml* or *config* file containing a dictionary called `deconv_args`, with the keys being the method name **in lowercase**.

```
# example.yaml
# To run: nextflow run main.nf -profile <PROFILE> -params-file example.yaml

deconv_args:
  stereoscope: "-stb 500 -scb 500"
  cell2location:
    build: "-t tech"
    fit: "-p 500"
  spatialdwls: "--n_topmarkers 50"
```

```
// example.config
// To run: nextflow run main.nf -profile <PROFILE> -c example.config

params.deconv_args = [stereoscope: "-stb 500 -scb 500",
                      cell2location: [build: "-t tech", fit: "-p 500"],
                      spatialdwls: "--n_topmarkers 50"]
```
Note that `params.sampleID` will automatically be fed to applicable methods once it is given as a pipeline parameter, so there is no need to feed it individually to methods, e.g., through `params.deconv_args.music`. We nonetheless included this information in the parameter list for the sake of completeness.



---

#### cell2location
Model building (`params.deconv_args.cell2location.build`)
- `-s`: metadata column in single-cell reference containing sample or batch information (default: None)
- `-t`: metadata column in single-cell reference containing multiplicative technical effects, such as platform effects (default: None)
- `-e`: number of epochs to train the model (default: 250)
- `-p`: number of samples to take from the posterior distribution (default: 1000)

Model fitting (`params.deconv_args.cell2location.fit`)
- `-n`: estimated number of cells per spot (default: 8)
- `-d`: within-experiment variation in RNA detection sensitivity (default: 200)
- `-e`: number of epochs to fit the model (default: 30000)
- `-p`: number of samples to take from the posterior distribution (default: 1000)

#### DestVI
Model building (`params.deconv_args.destvi.build`)
- `-e`: number of epochs to train the model (default: 300)
- `-n`: number of highly variable genes to use (default: 2000)

Model fitting (`params.deconv_args.destvi.fit`)
- `-e`: number of epochs to fit the model (default: 2500)
- `-b`: minibatch size to use during training (default: 128)

#### MuSiC
- `--sampleID`: metadata column in single-cell reference containing sample information (default: none)
- `--downsample_cells`: whether or not to downsample cells if the dense matrix is too large. The downsampling process will keep a maximum number of target_n_cells per cell type. (default: TRUE)
- `--target_n_cells`: the maximum number of cells per cell type after downsampling, if `downsample_cells=TRUE` (default: 10000)
- `--downsample_genes`: whether or not to downsample genes if the dense matrix is too large. The downsampling process will keep `n_hvgs` highly variable genes and filter out genes that are not expressed by at least `pct` fraction of cells per each cell type. (default: TRUE)
- `--n_hvgs`: number of highly variable genes to keep after gene downsampling (default: 3000)
- `--pct`: for the gene downsampling process, the fraction of cells per cell type in which genes have to be expressed in order to be kept (default: 0.1)
- `--assay_oi`: expression matrix to look for when downsampling genes (default: RNA)
- `--filter_spots`: minimum UMI count per spatial spot needed, otherwise it will be filtered out (default: none)

#### Seurat
- `--tech`: split the reference object based on this metadata column and integrate them, use only if the reference comes from different technologies (default: none)
- `--norm.method`: normalization method, either "vst" or "SCT" (default: vst)
- `--n_hvgs`: number of variable features to use (default: 2000)
- `--n_int_features`: number of variable features used for integration (default: 2000)
- `--npcs`: number of principal components (default: 30)
- `--dims`: number of dimensions for integration and other functions
- `--reduction`: passed to FindTransferAnchors, either "pcaproject, "lsiproject", "rpca", or "cca" (default: pcaproject)
- `--k.score`: passed to FindTransferAnchors (default: 30)
- `--k.weight`: passed to TransferData (default: 50)

#### spatialDWLS
- `--n_topmarkers`: number of top marker genes per cell type to use (default: 100)
- `--nn.dims`: number of PCs to use for nearest network calculation (default: 10)
- `--nn.k`: number of neighbors to use for nearest network calculation (default: 4)
- `--cluster.res`: Leiden cluster resolution (default: 0.4)
- `--cluster.n_iter`: number iterations during Leiden clustering (default: 1000)

#### SPOTlight
- `--logfc.threshold`: passed to Seurat::FindAllMarkers (default:0.25)
- `--min.pct`: passed to Seurat::FindAllMarkers (default: 0.1)
- `--conserve.memory`: passed to Seurat::SCTransform (default: FALSE)
- `--cl_n`: number of cells per cell type to use (default: 100)
- `--hvg`: number of HVGs to use (default: 3000)
- `--ntop`: how many of the marker genes to use, use all by default (default: NULL)
- `--transf`: whether to perform unit-variance per cell and spot (default: uv)
- `--method`: factorization method (default: nsNMF)
- `--min_cont`: remove cells contributing to a spot below a certain threshold (default: 0)

#### stereoscope
- `-scb`: batch size for single cell data set (default: None)
- `-sce`: number of epochs to be used in fitting of single cell data (default: 20000)
- `-stb`: batch size for spatial dataset (default: None)
- `-ste`: number of epochs to be used in fitting of spatial transcriptomics data. (default: 20000)
- `-kn`: keep noise (default: False)
- `-n`: only use top n mose highly expressed genes (default: None)
- `-fg`: filter Ribosomal Genes and MALAT1 (default: False)
- `-lr`: learning rate to be used (default: 0.01)
- `-mscc`: minimum number of counts for single cells to be included in the analysis (default: 0)
- `-mstc`: minimum number of counts for spots to be included in the analysis (default: 0)
- `-mc`: minimum number of cells for genes to be observed in the analysis (default: 0)
- `-ms`: minimum number of spots for genes to be observed in the analysis (default: 0)
- `-gl`: path to list of genes to use (default: None)
- `-hvg`: number of highly variable genes to use (default: None) 
- `-sub`: upper bound for single cell subsampling. (default: None)
- `-slb`: lower bound for single cell subsampling. (default: None)
- `-fb`: freeze beta parameter (default: False)

**Notes:** `-gl` and `-hvg` cannot be passed together, and `-hvg` will not work if there are spaces in your file name.

#### STRIDE
- `--gene-use`: location of the gene list file used to train the model, or "All" (longer runtime). By default the differential marker genes for each celltype will be used.
- `--markers`: number of differential marker genes to use (calculated using rank_genes_groups function)
- `--sc-scale-factor`: the scale factor for cell-level normalization (default: 75% quantile of nCount)
- `--st-scale-factor`: the scale factor for spot-level normalization (default: 75% quantile of nCount for ST)
- `--normalize`: flag indicating that the single-cell and spatial count matrices should be normalized (based on SD for each gene)
- `--ntopics`: number of topics to train and test the model. Multiple numbers should be separated by space (e.g., --ntopics 6 7 8). STRIDE will automatically select the optimal topic number (default: checks ntopics from number_of_celltypes to 3 $\times$ number_of_celltypes)

**Notes:** `--gene-use` and `--markers` cannot be passed together, and `--gene-use` will not work if there are spaces in your file name.


#### Tangram
- `-e`: number of epochs to train the model (default: 1000)
- `-n`: number of top marker genes to use (default: 100)
- `-m`: mode of running the deconvolution, either "cells", "clusters", or "constrained" (default: clusters)
- `-d`: density prior of the deconvolution, either "rna_count_based" or "uniform" (default: rna_count_based)
