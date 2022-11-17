import argparse as arp
import os
import sys

##### PARSING COMMAND LINE ARGUMENTS #####
prs = arp.ArgumentParser()

prs.add_argument('sc_data_path',
                    type = str, help = 'path to single cell h5ad count data')

prs.add_argument('sp_data_path', type = str, help = "path to a spatial dataset to get gene names" )

prs.add_argument('cuda_device', type = str, help = "index of cuda device ID or cpu")

prs.add_argument('-o','--out_dir', default = None,
                    type = str, help = 'directory for regression model')

prs.add_argument('-a','--annotation_column', default = 'celltype',
                type = str, help = 'column name for covariate')

prs.add_argument('-e', '--epochs', default=300, type = int,
                help = "number of epochs to train the model")

prs.add_argument('-n', '--n_hvgs', default=2000, type = int,
                help = "number of highly variable genes to use")

args = prs.parse_args()
cuda_device = args.cuda_device

assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device id"
assert os.path.exists(args.sc_data_path), "sc file not found"

if args.out_dir is None:
    output_folder = os.path.dirname(args.sc_data_path) + '/destvi_results/'
else:
    output_folder = args.out_dir

assert not os.path.exists(output_folder), "folder already exists"

print("Parameters\n==========")
print("Training epochs: {}".format(args.epochs))
print("==========")

##### MAIN CODE #####
if cuda_device.isdigit():
    os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device

import scanpy as sc
from scvi.model import CondSCVI, DestVI
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

## scRNA reference (raw counts)
print("Reading scRNA-seq data from " + args.sc_data_path + "...")
sc_adata = sc.read_h5ad(args.sc_data_path)

print("Reading in spatial data from " + args.sp_data_path + "...")
st_adata = sc.read_h5ad(args.sp_data_path)

if not all(sc_adata.var_names.isin(st_adata.var_names)):
    print("Subsetting single-cell data to match genes in spatial data...")
    print("Before subsetting: {} genes.".format(sc_adata.shape[1]))
    sc_adata = sc_adata[:, sc_adata.var_names.isin(st_adata.var_names)].copy()
    print("After subsetting: {} genes.".format(sc_adata.shape[1]))

# Filter genes
print("Before filtering: {} genes.".format(sc_adata.shape[1]))
G = args.n_hvgs
sc.pp.filter_genes(sc_adata, min_counts=10)
print("After filtering: {} genes.".format(sc_adata.shape[1]))
sc_adata.layers["counts"] = sc_adata.X.copy()

print("Subsetting on HVGs...")
sc.pp.highly_variable_genes(
    sc_adata,
    n_top_genes=G,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)
sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata

print("Single-cell data now has {} genes.".format(sc_adata.shape[1]))
print("Preparing anndata for the regression model...")
CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key=args.annotation_column)

# Run the model
sc_model = CondSCVI(sc_adata, weight_obs=False)
sc_model.view_anndata_setup()
sc_model.train(max_epochs=args.epochs)

# Save training figure
sc_model.history["elbo_train"].iloc[5:].plot()
plt.savefig("train.png")

try:
    sc_model.save(output_folder, save_anndata=True)
except ValueError:
    print("There seems to be an issue with the conversion. Renaming columns...")
    os.remove(output_folder + "/adata.h5ad")
    os.rmdir(output_folder)
    # sc_model.adata.var = sc_model.adata.var.drop(columns='_index')
    #adata_scrna_raw.__dict__['_raw'].__dict__['_var'] = adata_scrna_raw.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
        
    #sc_model.adata.__dict__['_raw'].__dict__['_var'] = sc_model.adata.__dict__['_raw'].__dict__['_var'].drop(columns='_index')
    sc_model.save(output_folder, save_anndata=True)


