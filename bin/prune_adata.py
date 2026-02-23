#!/usr/bin/env python3

import argparse
from runpy import run_path
import anndata as ad
import scanpy as sc
import numpy as np
import h5py

parser = argparse.ArgumentParser(description='Prune AnnData object.')
parser.add_argument('--raw_adata', type=str, help='Path to the input h5 file.')
parser.add_argument('--input_file', type=str, help='Path to the parameters file.')
parser.add_argument("--adata_mask", type=str, help="Path to the mask file.")
parser.add_argument('--n_pca_unintegrated', type=int, help='Number of PCA components to use for unintegrated data.')
args = parser.parse_args()

params = run_path(args.input_file)

with h5py.File(args.raw_adata) as file:
    obs = ad._io.specs.read_elem(file['obs'])
    
with h5py.File(args.raw_adata) as file:
    var = ad._io.specs.read_elem(file['var'])
    
adata = ad.AnnData(obs=obs, var=var)

if 'layer' in params['scvi_input']:
    with h5py.File(args.raw_adata) as file:
        adata.layers[params['scvi_input']['layer']] = ad._io.specs.read_elem(file[f"layers/{params['scvi_input']['layer']}"])
        adata.layers['nxf_norm'] = ad._io.specs.read_elem(file[f"layers/{params['scvi_input']['layer']}"])
else:
    with h5py.File(args.raw_adata) as file:
        adata.X = ad._io.specs.read_elem(file['X'])
        adata.layers['nxf_norm'] = ad._io.specs.read_elem(file['X'])

if args.adata_mask != '':
    adata = adata[:, adata.var[args.adata_mask]].copy()

if 'pre_integrated_embedding_obsm_key' not in params['scib_input']:
    sc.pp.normalize_total(adata, layer='nxf_norm')
    sc.pp.log1p(adata, layer='nxf_norm')
    sc.pp.pca(adata, n_comps=args.n_pca_unintegrated, layer='nxf_norm')

del adata.layers['nxf_norm']
del adata.uns
del adata.varm

adata.write_h5ad(f"pruned_adata_{args.adata_mask}.h5ad")

pca_adata = ad.AnnData(obs=adata.obs.copy(), var=adata.var.copy(), obsm={"X_pca": adata.obsm["X_pca"]})
pca_adata.write_h5ad(f"PCA_params_unintegrated_{args.adata_mask}.h5ad")
