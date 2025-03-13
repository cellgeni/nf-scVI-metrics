#!/usr/bin/env python3

import anndata as ad
import numpy as np
import scvi
import json
import pickle
import argparse
from runpy import run_path

parser = argparse.ArgumentParser(description="Run scVI model.")
parser.add_argument("--adata", type=str, help="Path to the input AnnData file.")
parser.add_argument("--input_file", type=str, help="Path to the arguments file.")
parser.add_argument("--params_file", type=str, help="Path to the model input JSON file.")
parser.add_argument("--check_val_every_n_epoch", type=int, help="Number of epochs to check validation loss.")
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

scvi_model.train(check_val_every_n_epoch = args.check_val_every_n_epoch, **params['train_input'])

np.save(f"scvi_{args.params_file}", scvi_model.get_latent_representation())

with open(f"history_{args.params_file}", "wb") as f:
    pickle.dump(scvi_model.history, f)
