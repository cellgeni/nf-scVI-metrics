#!/usr/bin/env python3

import sys
import h5py
import anndata as ad
import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from runpy import run_path
import argparse

parser = argparse.ArgumentParser(description="Run scIB metrics.")
parser.add_argument("--adata", type=str, help="Path to the input AnnData file.")
parser.add_argument("--input_file", type=str, help="Path to the arguments file.")
parser.add_argument("--scVI_embedding", type=str, help="Path to the embedding file.")
args = parser.parse_args()

adata = ad.read_h5ad(args.adata)
params = run_path(args.input_file)
temp = args.scVI_embedding

if 'layer' in params['param_input']:
        adata.X = adata.layers[params['param_input']['layer']]

embedding = np.load(temp)
# with h5py.File(temp) as file:
#     embedding = ad._io.specs.read_elem(file['obsm/scVI'])
temp_obsm = 'param_' + temp.split('_')[2].split('.')[0]
adata.obsm[temp_obsm] = embedding

scib_input = params['scib_input']
if 'pre_integrated_embedding_obsm_key' not in scib_input:
    scib_input['pre_integrated_embedding_obsm_key'] = 'X_pca'

scib_max_obs = 500000
if adata.shape[0] > scib_max_obs:
    sc.pp.subsample(adata, n_obs = scib_max_obs)

bm = Benchmarker(
    adata,
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    embedding_obsm_keys=[temp_obsm],
    n_jobs=2,
    **scib_input
)

bm.benchmark()

bm._results.to_csv(temp_obsm)
