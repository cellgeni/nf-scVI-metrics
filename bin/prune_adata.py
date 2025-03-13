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
parser.add_argument('--n_pca_unintegrated', type=int, help='Number of PCA components to use for unintegrated data.')
args = parser.parse_args()

params = run_path(args.input_file)

with h5py.File(args.raw_adata) as file:
    obs = ad._io.specs.read_elem(file['obs'])
    
with h5py.File(args.raw_adata) as file:
    var = ad._io.specs.read_elem(file['var'])
    
adata = ad.AnnData(obs=obs, var=var)

if 'layer' in params['param_input']:
    with h5py.File(args.raw_adata) as file:
        adata.layers[params['param_input']['layer']] = ad._io.specs.read_elem(file[f"layers/{params['param_input']['layer']}"])
        adata.layers['nxf_norm'] = ad._io.specs.read_elem(file[f"layers/{params['param_input']['layer']}"])
else:
    with h5py.File(args.raw_adata) as file:
        adata.X = ad._io.specs.read_elem(file['X'])
        adata.layers['nxf_norm'] = ad._io.specs.read_elem(file['X'])

if 'pre_integrated_embedding_obsm_key' not in params['scib_input']:
    sc.pp.normalize_total(adata, layer='nxf_norm')
    sc.pp.log1p(adata, layer='nxf_norm')
    sc.pp.pca(adata, n_comps=args.n_pca_unintegrated, layer='nxf_norm')

del adata.layers['nxf_norm']
del adata.uns
del adata.varm

adata.write_h5ad("pruned_adata")

np.save(f"PCA_params_unintegrated", adata.obsm['X_pca'])
