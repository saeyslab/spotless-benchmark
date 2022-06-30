
import argparse as arp
import os
import sys

##### PARSING COMMAND LINE ARGUMENTS #####
prs = arp.ArgumentParser()

prs.add_argument('sc_data_path',
                    type = str, help = 'path to single cell h5ad count data')

prs.add_argument('sp_data_path', type = str, help = 'path to h5ad spatial data')

prs.add_argument('cuda_device', type = str, help = "index of cuda device ID or cpu")

prs.add_argument('-o','--output_file', default = "proportions.tsv",
                    type = str, help = 'proportions output file')

prs.add_argument('-a','--annotation_column', default = 'celltype',
                type = str, help = 'column name for covariate')

prs.add_argument('-c', '--cell_count_file', type = str, default = None,
                help = "path to file with cell count information")

prs.add_argument('-e', '--epochs', default=1000, type = int,
                help = "number of epochs to train the model")

prs.add_argument('-n', '--n_top_markers', default=100, type = int,
                help = "number of top marker genes to use")

prs.add_argument('-m', '--mode', default="clusters", type = str,
                help = "mode of running the deconvolution (cells, clusters, or constrained)")

prs.add_argument('-d', '--density_prior', default="rna_count_based", type = str,
                help = "density prior of the deconvolution")

args = prs.parse_args()
cuda_device = args.cuda_device

assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device id"
assert os.path.exists(args.sc_data_path), "sc file not found"
assert os.path.exists(args.sp_data_path), "sp file not found"
if args.cell_count_file:
    assert os.path.exists(args.cell_count_file), "cell count file not found"
assert args.mode in ['cells', 'clusters', 'constrained'], "deconvolution mode must be either, cells, clusters, or constrained"
assert args.density_prior in ['rna_count_based', 'uniform'], "density prior must be either rna_count_based or uniform"

print("Parameters\n==========")
print("Training epochs: {}".format(args.epochs))
print("Deconvolution mode: {}".format(args.mode))
print("Density prior: {}".format(args.density_prior))
print("Number of top markers: {}".format(args.n_top_markers))
print("==========")

##### MAIN CODE #####
import scanpy as sc
import tangram as tg
import pandas as pd
import numpy as np

if cuda_device.isdigit():
    os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device

## scRNA reference (raw counts)
print("Reading scRNA-seq data from " + args.sc_data_path + "...")
sc_adata = sc.read_h5ad(args.sc_data_path)

print("Reading spatial data from " + args.sp_data_path + "...")
sp_adata = sc.read_h5ad(args.sp_data_path)

print("Normalizing and logarithmizing data...")
sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)

print("Getting marker genes...")
sc.tl.rank_genes_groups(sc_adata, groupby=args.annotation_column, use_raw=False)
markers_df = pd.DataFrame(sc_adata.uns["rank_genes_groups"]["names"]).iloc[0:args.n_top_markers, :]
markers = list(np.unique(markers_df.melt().value.values))
tg.pp_adatas(sc_adata, sp_adata, genes=markers)

# Initialize hyperparameters
deconv_mode = args.mode
density_prior = args.density_prior
target_counts = None

# In this mode, the number of single-cell profiles to be mapped is constrained
# We either get the cell counts from the ground truth or from segmenting the H&E image
if deconv_mode == "constrained":
    cell_counts = pd.read_csv(args.cell_count_file)

    # If no ground truth or image slot was found in getCellComposition.R
    if cell_counts.iloc[0][0] == 'no_counts_found':
        print("Changing 'mode' to 'clusters' and density prior to {}...".format(density_prior))
        deconv_mode = "clusters"

    else:
        # If an image slot was found, we use squidpy to segment the image
        if cell_counts.iloc[0][0] == 'segment':
            import squidpy as sq
            print("Cell counts will be obtained from segmenting the h5ad file...")
            print("Creating squidpy ImageContainer from hires image...")
            img_data = sp_adata.uns['spatial'][list(sp_adata.uns['spatial'].keys())[0]]
            img = sq.im.ImageContainer(img_data['images']['hires'],
                                    layer = "image",
                                    scale = img_data['scalefactors']['tissue_hires_scalef'])
            
            print("Preprocessing and segmenting image...")
            sq.im.process(img=img, layer="image", method="smooth")
            sq.im.segment(img=img, layer="image_smooth",
                        method="watershed", channel=0)

            # Define image layer to use for segmentation
            features_kwargs = {"segmentation": { "label_layer": "segmented_watershed",
                                                "props": ["label", "centroid"],
                                                "channels": [1, 2] }}
            print ("Calculating segmentation features...")
            sq.im.calculate_image_features(
                sp_adata, img,
                layer="image",
                key_added="image_features",
                features_kwargs=features_kwargs,
                features="segmentation",
                mask_circle=True
            )

            sp_adata.obs["cell_count"] = sp_adata.obsm["image_features"]["segmentation_label"]
        
        # Otherwise, the cell counts exist as a ground truth
        else:
            print("Cell counts was computed in previous step...")
            sp_adata.obs["cell_count"] = cell_counts
        
        # Calculate total cells and cell density per spot
        target_counts = sp_adata.obs.cell_count.sum()
        density_prior = np.array(sp_adata.obs.cell_count) / sp_adata.obs.cell_count.sum()

print("Final Parameters\n==========")
print("Deconvolution mode: {}".format(deconv_mode))
print("Density prior: {}".format(density_prior if type(density_prior) == str else "custom array"))
print("Target counts: {}".format(target_counts))
print("==========")


# Deconvolution
ad_map = tg.map_cells_to_space(
    sc_adata, sp_adata,
    mode=deconv_mode,
    target_count=target_counts,
    density_prior=density_prior,
    cluster_label=None if deconv_mode != "clusters" else args.annotation_column,
    num_epochs=args.epochs,
    device='cuda:'+cuda_device if cuda_device.isdigit() else cuda_device
)

tg.project_cell_annotations(ad_map, sp_adata, annotation=args.annotation_column)

deconv_matrix = sp_adata.obsm['tangram_ct_pred']/np.sum(sp_adata.obsm['tangram_ct_pred'], axis=1)[:, None]
deconv_matrix.to_csv(args.output_file, sep="\t")
