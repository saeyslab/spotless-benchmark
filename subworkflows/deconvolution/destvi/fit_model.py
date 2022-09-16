#!/usr/bin/env python3
import argparse as arp
import os
import sys

##### PARSING COMMAND LINE ARGUMENTS #####
prs = arp.ArgumentParser()

prs.add_argument('sp_data_path', type = str, help = 'path of spatial data')

prs.add_argument('cuda_device', type = str, help = "index of cuda device ID or cpu")

prs.add_argument('-m', '--model_path', default = "model/", type = str, help = 'path to regression model')

prs.add_argument('-o','--out_dir', default = os.getcwd() ,
                    type = str, help = 'model and proportions output directory')

prs.add_argument('-e', '--epochs', default=2500, type = int, help = "number of epochs to fit the model")

prs.add_argument('-b', '--batch_size', default=128, type=int, help = "minibatch size to use during training.")

args = prs.parse_args()

cuda_device = args.cuda_device
sp_data_path = args.sp_data_path
output_folder = args.out_dir
    
assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device input"

print("Parameters\n==========")
print("Training epochs: {}".format(args.epochs))
print("Batch size: {}".format(args.batch_size))
print("==========")

##### MAIN PART #####
if cuda_device.isdigit():
    os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device

import scanpy as sc
import numpy as np
from scvi.model import CondSCVI, DestVI
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

print("Reading in spatial data from " + sp_data_path + "...")
st_adata = sc.read_h5ad(sp_data_path)
st_adata.layers["counts"] = st_adata.X.copy()

print("Reading in the sc model...")
sc_model = CondSCVI.load(args.model_path)

if st_adata.shape[1] != sc_model.adata.shape[1]:
    print("The number of genes do not match. Subsetting spatial data...")
    st_adata = st_adata[:, st_adata.var_names.isin(sc_model.adata.var_names)].copy()

# Prepare anndata
print("Setting up spatial model...")
DestVI.setup_anndata(st_adata, layer="counts")

# Set up model
st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.view_anndata_setup()
st_model.train(max_epochs=args.epochs, batch_size=args.batch_size)

# Save training figure
st_model.history["elbo_train"].iloc[5:].plot()
plt.savefig("train.png")

# Export proportion file
st_model.get_proportions().to_csv(os.path.join(output_folder, 'proportions.tsv'), sep="\t")


    
