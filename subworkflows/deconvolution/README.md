## Method parameters
Parameters for each method can be provided in a *yaml* or *config* file containing a dictionary called `deconv_args`, with the keys being the method name in lowercase.
Note that the `params.sampleID` will automatically be fed to applicable methods once it is given as a pipeline parameter, so there is no need to feed it individually to methods, e.g., within
`params.deconv_args.music`. We will however mention below which method can handle sample information for the sake of completeness.

```
# example.yaml
# To run: nextflow run main.nf -profile <PROFILE> -params-file example.yaml

deconv_args:
  stereoscope: "-stb 500 -scb 500"
  cell2location:
    build: "-t tech"
    fit: "-p 500"
```

```
// example.config
// To run: nextflow run main.nf -profile <PROFILE> -c example.config

params.deconv_args = [stereoscope: "-stb 500 -scb 500",
                      cell2location: [build: "-t tech", fit: "-p 500"]]
```

#### cell2location

#### DestVI

#### DSTG

#### MuSiC
- `--sampleID`: metadata column in the single-cell reference file containing sample information (default: none)
- `downsample_cells`: whether or not to downsample cells if the dense matrix is too large. The downsampling process will keep a maximum number of target_n_cells per cell type. (default: TRUE)
- `target_n_cells`: the maximum number of cells per cell type after downsampling, if `downsample_cells=TRUE` (default: 10000)
- `downsample_genes`: whether or not to downsample genes if the dense matrix is too large. The downsampling process will keep `n_hvgs` highly variable genes and filter out genes that are not expressed by at least `pct` fraction of cells per each cell type. (default: TRUE)
- `n_hvgs`: number of highly variable genes to keep after gene downsampling (default: 3000)
- `pct`: for the gene downsampling process, the fraction of cells per cell type in which genes have to be expressed in order to be kept (default: 0.1)
- `assay_oi`: expression matrix to look for when downsampling genes (default: "RNA")
- `filter_spots`: minimum UMI count per spatial spot needed, otherwise it will be filtered out (default: "none")

#### RCTD
No adjustable parameters

#### spatialDWLS

#### SPOTlight

#### stereoscope
