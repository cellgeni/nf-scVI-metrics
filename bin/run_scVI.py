import sys
import anndata as ad
import scanpy as sc
import scvi
import json
import pickle

from runpy import run_path

adata_path = sys.argv[1]
args = run_path(sys.argv[2])

with open(sys.argv[3]) as f:
    model_input = json.load(f)

adata = ad.read_h5ad(adata_path)

scvi.model.SCVI.setup_anndata(
    adata,
    **args['param_input']
)

scvi_model = scvi.model.SCVI(adata,
                             **model_input)

scvi_model.train(check_val_every_n_epoch = 5, **args['train_input'])

import numpy as np
np.save(f"param_{sys.argv[3]}", scvi_model.get_latent_representation())

# out_adata = ad.AnnData(obs = adata.obs, var = adata.var)
# out_adata.obsm['scVI'] = scvi_model.get_latent_representation()
# out_adata.write_h5ad(f"adata_{sys.argv[3]}")

with open(f"history_{sys.argv[3]}", "wb") as f:
    pickle.dump(scvi_model.history, f)
