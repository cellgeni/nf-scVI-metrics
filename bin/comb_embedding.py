import sys
import numpy as np
import anndata as ad

adata = ad.read_h5ad(sys.argv[1])

for i in sorted(sys.argv[2].split()):
    adata.obsm[f"scvi_{i.split('_')[3]}"] = np.load(i)

if sys.argv[3]:
    for i in sorted(sys.argv[4].split()):
        adata.obsm[f"umap_{i.split('_')[3]}"] = np.load(i)

adata.write_h5ad("adata_embedding.h5ad")