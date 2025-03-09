import sys
from runpy import run_path
import anndata as ad
import scanpy as sc
import numpy as np
import h5py

args = run_path(sys.argv[2])

with h5py.File(sys.argv[1]) as file:
    obs = ad._io.specs.read_elem(file['obs'])
    
with h5py.File(sys.argv[1]) as file:
    var = ad._io.specs.read_elem(file['var'])
    
adata = ad.AnnData(obs = obs, var = var)

if 'layer' in args['param_input']:
    with h5py.File(sys.argv[1]) as file:
        adata.layers[args['param_input']['layer']] = ad._io.specs.read_elem(file[f"layers/{args['param_input']['layer']}"])
        adata.layers['nxf_norm'] = ad._io.specs.read_elem(file[f"layers/{args['param_input']['layer']}"])
else:
    
    with h5py.File(sys.argv[1]) as file:
        adata.X = ad._io.specs.read_elem(file['X'])
        adata.layers['nxf_norm'] = ad._io.specs.read_elem(file['X'])

if 'pre_integrated_embedding_obsm_key' not in args['scib_input']:
    sc.pp.normalize_total(adata, layer = 'nxf_norm')
    sc.pp.log1p(adata, layer = 'nxf_norm')
    sc.pp.pca(adata, n_comps = 30, layer = 'nxf_norm')

del adata.layers['nxf_norm']
del adata.uns
del adata.varm

adata.write_h5ad("input_adata")

np.save(f"PCA_prune_data_unintegrated", adata.obsm['X_pca'])
