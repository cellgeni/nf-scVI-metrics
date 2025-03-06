import sys
from runpy import run_path
import anndata as ad
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
else:
    
    with h5py.File(sys.argv[1]) as file:
        adata.X = ad._io.specs.read_elem(file['X'])

adata.write_h5ad("input_adata")
