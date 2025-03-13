#!/usr/bin/env python3

import argparse
import anndata as ad
import scanpy as sc
import numpy as np

parser = argparse.ArgumentParser(description='Run UMAP on scVI embeddings.')
parser.add_argument('--scVI_embedding', type=str, help='Path to the input numpy file.')
parser.add_argument('--seed', type=int, help='Random seed for UMAP.')
args = parser.parse_args()

adata = ad.AnnData(np.load(args.scVI_embedding))
sc.pp.neighbors(adata, use_rep='X')
sc.tl.umap(adata, random_state=args.seed)
np.save(f"umap_{args.scVI_embedding}", adata.obsm['X_umap'])
