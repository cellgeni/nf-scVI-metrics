import sys
import anndata as ad
import scanpy as sc

adata_name = sys.argv[1]

adata = ad.read_h5ad(adata_name)

sc.pp.neighbors(adata, use_rep='scVI')
sc.tl.umap(adata, random_state=123)

adata.write_h5ad(f"umap_{adata_name}")
