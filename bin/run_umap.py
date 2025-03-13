#!/usr/bin/env python3

import sys
import argparse
import anndata as ad
import scanpy as sc
import numpy as np

parser = argparse.ArgumentParser(description='Run UMAP on scVI embeddings.')
parser.add_argument('--scVI_embedding', type=str, help='Path to the input numpy file.')
args = parser.parse_args()

adata = ad.AnnData(np.load(args.scVI_embedding))
sc.pp.neighbors(adata, use_rep='X')
sc.tl.umap(adata, random_state=123)
np.save(f"umap_{args.scVI_embedding}", adata.obsm['X_umap'])
