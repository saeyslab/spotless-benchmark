import scanpy as sc
import sys
import re

args = " ".join(sys.argv[1:])

m = re.search('(?<=-hvg )[0-9]+', args)
n_hvgs = m.group(0)

# This will not work if you have spaces in your input file
m = re.search('(?<=--sc_cnt )[^ ]+', args)
sc_path = m.group(0)

print("Computing {} HVGs from {}...".format(n_hvgs, sc_path))

adata = sc.read_h5ad(sc_path)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=int(n_hvgs), subset=True)

adata.var_names.to_frame().to_csv("hvgs.txt", header=False, index=False)