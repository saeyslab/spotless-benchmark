#!/usr/bin/env python3
import argparse as arp
import os

def main():
    ##### PARSING COMMAND LINE ARGUMENTS #####
    prs = arp.ArgumentParser()
    
    prs.add_argument('sp_data_path', type = str, help = 'path of spatial data')

    prs.add_argument('model_path', type = str, help = 'path to regression model')
    
    prs.add_argument('cuda_device', type = str, help = "index of cuda device ID or cpu")

    prs.add_argument('-o','--out_dir', default = os.getcwd() ,
                     type = str, help = 'model and proportions output directory')

    prs.add_argument('-e', '--epochs', default=30000, type = int, help = "number of epochs to fit the model")

    prs.add_argument('-p', '--posterior_sampling', default=1000, type = int, help = "number of samples to take from the posterior distribution")

    prs.add_argument('-n', '--n_cells_per_location', default=8, type = int, help = "estimated number of cells per spot")

    prs.add_argument('-d', '--detection_alpha', default=200, type = int, help = "within-experiment variation in RNA detection sensitivity")
    
    args = prs.parse_args()
    
    cuda_device = args.cuda_device
    sp_data_path = args.sp_data_path
    output_folder = args.out_dir
        
    assert (cuda_device.isdigit() or cuda_device == "cpu"), "invalid device input"
    
    print("Parameters\n==========")
    print("Detection alpha: {}\nCells per location: {}".format(args.detection_alpha, args.n_cells_per_location))
    print("Training epochs: {}\nPosterior sampling: {}".format(args.epochs, args.posterior_sampling))
    print("==========")

    ##### MAIN PART #####
    if cuda_device.isdigit():
        os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device

    import sys
    import scanpy as sc
    import anndata
    import pandas as pd
    import numpy as np

    import cell2location
    import scvi

    import matplotlib as mpl
    from matplotlib import rcParams
    rcParams['pdf.fonttype'] = 42
    import matplotlib.pyplot as plt
    import seaborn as sns

    # silence scanpy that prints a lot of warnings
    import warnings
    warnings.filterwarnings('ignore')

    print("Reading in spatial data from " + sp_data_path + "...")
    adata = sc.read_h5ad(sp_data_path)
    adata.var['SYMBOL'] = adata.var_names

    # mitochondria-encoded (MT) genes should be removed for spatial mapping
    adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata = adata[:, ~adata.var['mt'].values]

    adata_vis = adata.copy()
    adata_vis.raw = adata_vis

    print("Reading in the model...")
    adata_scrna_raw = sc.read(args.model_path)
    
    # Export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_scrna_raw.varm.keys():
        inf_aver = adata_scrna_raw.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                        for i in adata_scrna_raw.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_scrna_raw.var[[f'means_per_cluster_mu_fg_{i}' 
                                        for i in adata_scrna_raw.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_scrna_raw.uns['mod']['factor_names']

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    scvi.data.setup_anndata(adata=adata_vis)

    # Create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver, 
        # the expected average cell abundance: tissue-dependent 
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=args.n_cells_per_location,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection (using default here):
        detection_alpha=args.detection_alpha
    ) 

    mod.train(max_epochs=args.epochs, 
            # train using full data (batch_size=None)
            batch_size=None, 
            # use all data points in training because 
            # we need to estimate cell abundance at all locations
            train_size=1,
            use_gpu=cuda_device.isdigit())

    # Export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': args.posterior_sampling,
        'batch_size': mod.adata.n_obs, 'use_gpu': cuda_device.isdigit()}
    )

    # Save model and anndata object with results
    mod.save(output_folder , overwrite=True)
    adata_vis.write(os.path.join(output_folder, 'sp.h5ad'))

    # Export proportion file, but first rename columns and divide by rowSums
    props = adata_vis.obsm['q05_cell_abundance_w_sf']
    props = props.rename(columns={x:x.replace("q05cell_abundance_w_sf_", "") for x in props.columns})
    props = props.div(props.sum(axis=1), axis='index')
    props.to_csv(os.path.join(output_folder, 'proportions.tsv'), sep="\t")

    # df = pd.DataFrame(data=np.random.normal(size=(10,10)),
    #                     index=["row"+str(i) for i in range(10)],
    #                     columns=["col"+str(i) for i in range(10)])
    # df.to_csv(os.path.join(output_folder, 'proportions.tsv'), sep="\t")

        
if __name__ == '__main__':
    main()