#!/usr/bin/env python3

import argparse
import anndata as ad
import scanpy as sc
import numpy as np
from pathlib import Path

parser = argparse.ArgumentParser(description='Run UMAP on scVI embeddings.')
parser.add_argument('--embedding_h5ad', type=str, help='Path to the embedding h5ad file.')
parser.add_argument('--seed', type=int, help='Random seed for UMAP.')
args = parser.parse_args()

emb = ad.read_h5ad(args.embedding_h5ad)
if "X_scVI" in emb.obsm:
    use_rep = "X_scVI"
elif "X_scANVI" in emb.obsm:
    use_rep = "X_scANVI"
else:
    raise KeyError("No supported embedding key found in embedding h5ad (expected X_scVI or X_scANVI).")

adata = ad.AnnData(emb.obsm[use_rep])
sc.pp.neighbors(adata, use_rep='X')
sc.tl.umap(adata, random_state=args.seed)
emb.obsm["X_umap"] = adata.obsm["X_umap"]
out_name = f"umap_{Path(args.embedding_h5ad).stem}.h5ad"
emb.write_h5ad(out_name)
