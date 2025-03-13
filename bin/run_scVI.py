#!/usr/bin/env python3

import sys
import anndata as ad
import scanpy as sc
import scvi
import json
import pickle
import argparse
from runpy import run_path

parser = argparse.ArgumentParser(description="Run scVI model.")
parser.add_argument("--adata", type=str, help="Path to the input AnnData file.")
parser.add_argument("--input_file", type=str, help="Path to the arguments file.")
parser.add_argument("--params_file", type=str, help="Path to the model input JSON file.")
args = parser.parse_args()

adata_path = args.adata
params = run_path(args.input_file)

with open(args.params_file) as f:
    model_input = json.load(f)

adata = ad.read_h5ad(adata_path)

scvi.model.SCVI.setup_anndata(
    adata,
    **params['param_input']
)

scvi_model = scvi.model.SCVI(adata,
                             **model_input)

scvi_model.train(check_val_every_n_epoch = 5, **params['train_input'])

import numpy as np
np.save(f"scvi_{args.params_file}", scvi_model.get_latent_representation())

# out_adata = ad.AnnData(obs = adata.obs, var = adata.var)
# out_adata.obsm['scVI'] = scvi_model.get_latent_representation()
# out_adata.write_h5ad(f"adata_{args.params_file}")

with open(f"history_{args.params_file}", "wb") as f:
    pickle.dump(scvi_model.history, f)
