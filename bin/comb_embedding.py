#!/usr/bin/env python3

import numpy as np
import anndata as ad
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Combine embeddings into an AnnData object.")
parser.add_argument('--raw_adata', type=str, help="Path to the input AnnData file.")
parser.add_argument('--embeddings', type=str, help="Space-separated list of embedding files.")
args = parser.parse_args()

adata = ad.read_h5ad(args.raw_adata)

for emb_path in sorted(args.embeddings.split()):
    path = Path(emb_path)
    if path.suffix == ".h5ad":
        emb = ad.read_h5ad(path)
        if "X_scVI" in emb.obsm:
            adata.obsm["X_scVI"] = emb.obsm["X_scVI"]
        if "X_scANVI" in emb.obsm:
            adata.obsm["X_scANVI"] = emb.obsm["X_scANVI"]
        if "X_pca" in emb.obsm:
            adata.obsm["X_pca"] = emb.obsm["X_pca"]
        if "X_umap" in emb.obsm:
            adata.obsm["X_umap"] = emb.obsm["X_umap"]
    else:
        adata.obsm[path.stem] = np.load(path)

adata.write_h5ad("adata_embedding.h5ad")
