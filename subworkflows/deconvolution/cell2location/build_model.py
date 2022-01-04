#!/usr/bin/env python3

import argparse as arp
import sys
import os
import matplotlib as mpl
mpl.use('Agg')

def main():
    ##### PARSING COMMAND LINE ARGUMENTS #####
    prs = arp.ArgumentParser()
    
    prs.add_argument('sc_data_path',
                     type = str, help = 'path to single cell h5ad count data')
    
    prs.add_argument('cuda_device', type = str, help = "index of cuda device ID or cpu")

    prs.add_argument('-o','--out_dir', default = None,
                     type = str, help = 'directory for regression model')

    prs.add_argument('-a','--annotation_column', default = 'celltype',
                 type = str, help = 'column name for covariate')

    prs.add_argument('-s', '--sample_column', default = None,
                 type = str, help = 'column containing sample id or batch information')
    
    prs.add_argument('-t', '--tech_column', default = None, nargs='+',
                 type = str, help = "multiplicative technical effects, such as platform effects")
    
    prs.add_argument('-e', '--epochs', default=250, type = int, help = "number of epochs to train the model")

    prs.add_argument('-p', '--posterior_sampling', default=1000, type = int, help = "number of samples to take from the posterior distribution")

    
    args = prs.parse_args()
    cuda_device = args.cuda_device

    assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device id"
    assert os.path.exists(args.sc_data_path), "sc file not found"

    if args.out_dir is None:
        output_folder = os.path.dirname(args.sc_data_path) + '/cell2location_results/'
    else:
        output_folder = args.out_dir

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print("Parameters\n==========")
    print("Training epochs: {}\nPosterior sampling: {}".format(args.epochs, args.posterior_sampling))
    print("==========")
    
    ##### MAIN CODE #####
    if cuda_device.isdigit():
        os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device

    import anndata
    import scanpy as sc
    import pandas as pd
    import numpy as np

    import cell2location
    import scvi
    
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams['pdf.fonttype'] = 42
    import seaborn as sns

    # silence scanpy that prints a lot of warnings
    import warnings
    warnings.filterwarnings('ignore')

    ## scRNA reference (raw counts)
    print("Reading scRNA-seq data from " + args.sc_data_path + "...")
    adata_scrna_raw = anndata.read_h5ad(args.sc_data_path)
    adata_scrna_raw.var['SYMBOL'] = adata_scrna_raw.var_names

    # Filter genes
    print("Before filtering: {} cells and {} genes.".format(*adata_scrna_raw.shape))
    from cell2location.utils.filtering import filter_genes
    selected = filter_genes(adata_scrna_raw, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    adata_scrna_raw = adata_scrna_raw[:, selected].copy()
    print("After filtering: {} cells and {} genes.".format(*adata_scrna_raw.shape))

    print("Preparing anndata for the regression model...")
    scvi.data.setup_anndata(adata=adata_scrna_raw, 
                        # 10X reaction / sample / batch
                        batch_key=args.sample_column, 
                        # cell type, covariate used for constructing signatures
                        labels_key=args.annotation_column, 
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=args.tech_column)
    
    # Run the model - use all data for training (validation not implemented yet, train_size=1)
    from cell2location.models import RegressionModel
    mod = RegressionModel(adata_scrna_raw) 
    mod.train(max_epochs=args.epochs, batch_size=2500, train_size=1, lr=0.002, use_gpu=cuda_device.isdigit())

    # Export the estimated cell abundance (summary of the posterior distribution).
    adata_scrna_raw = mod.export_posterior(
        adata_scrna_raw, sample_kwargs={'num_samples': args.posterior_sampling, 'batch_size': 2500, 'use_gpu': cuda_device.isdigit()}
    )

    # Save model and anndata object with results
    mod.save(output_folder, overwrite=True)

    try:
        adata_scrna_raw.write(output_folder + '/sc.h5ad')
    except ValueError:
        print("There seems to be an issue with the conversion. Renaming columns...")
        os.remove(output_folder + '/sc.h5ad')
        adata_scrna_raw.__dict__['_raw'].__dict__['_var'] = adata_scrna_raw.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
        adata_scrna_raw.write(output_folder + '/sc.h5ad')

    # with open(output_folder + '/sc.h5ad', "w") as f:
    #    print("hello world", file=f)


if __name__ == '__main__':
    main()