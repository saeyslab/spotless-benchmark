import scanpy as sc
import sys


filename = sys.argv[1]
annot = sys.argv[2]

adata = sc.read_h5ad(filename)
print("Writing the celltype annotation as annot.txt...")
adata.obs[annot].to_csv('annot.txt', sep='\t', header=False)

if len(sys.argv) > 3:
    print("Saving dummy ST file containing first 500 rows of sc file...")
    dummy_st = adata[:500,:]     
    try:
        dummy_st.write_h5ad("dummy_st.h5ad")
    except ValueError:
        print("There seems to be an issue with writing the file. Renaming columns...")
        dummy_st.__dict__['_raw'].__dict__['_var'] = dummy_st.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
        dummy_st.write_h5ad('dummy_st.h5ad')
