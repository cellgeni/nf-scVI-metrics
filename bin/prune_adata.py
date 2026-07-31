#!/usr/bin/env python3

import argparse
import anndata as ad
import scanpy as sc
import numpy as np
import h5py
import yaml

parser = argparse.ArgumentParser(description='Prune AnnData object.')
parser.add_argument('--raw_adata', type=str, help='Path to the input h5 file.')
parser.add_argument('--input_file', type=str, help='Path to the parameters file.')
parser.add_argument("--adata_mask", type=str, help="Path to the mask file.")
parser.add_argument('--n_pca_unintegrated', type=int, help='Number of PCA components to use for unintegrated data.')
args = parser.parse_args()

with open(args.input_file, "r", encoding="utf-8") as handle:
    params = yaml.safe_load(handle)
preprocess_config = params["preprocess"]
scib_config = params["scib"]["metrics"]

with h5py.File(args.raw_adata) as file:
    obs = ad._io.specs.read_elem(file['obs'])
    
with h5py.File(args.raw_adata) as file:
    var = ad._io.specs.read_elem(file['var'])
    
adata = ad.AnnData(obs=obs, var=var)

if "layer" in preprocess_config:
    with h5py.File(args.raw_adata) as file:
        adata.layers[preprocess_config["layer"]] = ad._io.specs.read_elem(file[f"layers/{preprocess_config['layer']}"])
        adata.layers['nxf_norm'] = ad._io.specs.read_elem(file[f"layers/{preprocess_config['layer']}"])
else:
    with h5py.File(args.raw_adata) as file:
        adata.X = ad._io.specs.read_elem(file['X'])
        adata.layers['nxf_norm'] = ad._io.specs.read_elem(file['X'])

if args.adata_mask != '':
    adata = adata[:, adata.var[args.adata_mask]].copy()

if 'pre_integrated_embedding_obsm_key' not in scib_config:
    sc.pp.normalize_total(adata, layer='nxf_norm')
    sc.pp.log1p(adata, layer='nxf_norm')
    sc.pp.pca(adata, n_comps=args.n_pca_unintegrated, layer='nxf_norm')

del adata.layers['nxf_norm']
del adata.uns
del adata.varm

adata.write_h5ad(f"pruned_adata_{args.adata_mask}.h5ad")


np.save(f"PCA_params_unintegrated_{args.adata_mask}", adata.obsm['X_pca'])