#!/usr/bin/env python3

import sys
import numpy as np
import anndata as ad
import argparse

parser = argparse.ArgumentParser(description="Combine embeddings into an AnnData object.")
parser.add_argument('--raw_adata', type=str, help="Path to the input AnnData file.")
parser.add_argument('--embeddings', type=str, help="Space-separated list of embedding files.")
args = parser.parse_args()

adata = ad.read_h5ad(args.raw_adata)

for i in sorted(args.embeddings.split()):
    adata.obsm[i.split('.')[0]] = np.load(i)

adata.write_h5ad("adata_embedding.h5ad")