import scanpy as sc
import pandas as pd
import sys
import re

args = " ".join(sys.argv[1:])

m = re.search('(?<=--markers )[0-9]+', args)
n_markers = int(m.group(0))

# This will not work if you have spaces in your input file
m = re.search('(?<=--sc-count )[^ ]+', args)
sc_path = m.group(0)

m = re.search('(?<=--annot )[^ ]+', args)
celltype = m.group(0)

print("Computing {} markers from {}...".format(n_markers, sc_path))

adata = sc.read_h5ad(sc_path)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.tl.rank_genes_groups(adata, celltype, method='wilcoxon', n_genes=n_markers)

markers = []
for i in adata.uns['rank_genes_groups']['names']:
    markers.extend(i)
    
pd.DataFrame(set(markers)).to_csv("markers.txt", header=False, index=False)