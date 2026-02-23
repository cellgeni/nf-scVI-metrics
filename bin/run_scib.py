#!/usr/bin/env python3

import anndata as ad
import scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from runpy import run_path
import argparse

parser = argparse.ArgumentParser(description="Run scIB metrics.")
parser.add_argument("--adata", type=str, help="Path to the input AnnData file.")
parser.add_argument("--input_file", type=str, help="Path to the arguments file.")
parser.add_argument("--embedding_h5ad", type=str, help="Path to the embedding h5ad file.")
parser.add_argument("--scib_max_obs", type=int, help="Maximum number of observations to use in scIB.")
parser.add_argument("--n_cpu", type=int, help="Number of CPUs to use in scIB.")
parser.add_argument("--scib_label_key", type=str, default="", help="Override label_key for scib_input.")
args = parser.parse_args()

adata = ad.read_h5ad(args.adata)
params = run_path(args.input_file)
embedding_path = args.embedding_h5ad

if 'layer' in params['scvi_input']:
        adata.X = adata.layers[params['scvi_input']['layer']]

def infer_obsm_key(emb_adata) -> str:
    if "X_scVI" in emb_adata.obsm:
        return "X_scVI"
    if "X_scANVI" in emb_adata.obsm:
        return "X_scANVI"
    if "X_pca" in emb_adata.obsm:
        return "X_pca"
    raise KeyError("No supported embedding key found in embedding h5ad (expected X_scVI, X_scANVI, or X_pca).")

emb_adata = ad.read_h5ad(embedding_path)
temp_obsm = infer_obsm_key(emb_adata)
temp_name = embedding_path.replace('scanvi_model_', '').replace(".h5ad", "")
adata.obsm[temp_name] = emb_adata.obsm[temp_obsm]

if not adata.obs_names.equals(emb_adata.obs_names):
    raise ValueError("Embedding h5ad obs_names do not match adata obs_names.")
for col in emb_adata.obs.columns:
    adata.obs[col] = emb_adata.obs[col]

scib_input = params['scib_input']
if args.scib_label_key:
    scib_input['label_key'] = args.scib_label_key
if 'pre_integrated_embedding_obsm_key' not in scib_input:
    scib_input['pre_integrated_embedding_obsm_key'] = 'X_pca'

label_key = scib_input.get('label_key')
if label_key and label_key not in emb_adata.obs.columns:
    raise KeyError(f"label_key '{label_key}' not found in embedding h5ad obs.")

scib_max_obs = args.scib_max_obs
if adata.shape[0] > scib_max_obs:
    sc.pp.subsample(adata, n_obs = scib_max_obs)

bm = Benchmarker(
    adata,
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    embedding_obsm_keys=[temp_name],
    n_jobs=args.n_cpu,
    **scib_input
)

bm.benchmark()

bm._results.to_csv('X_' + temp_name + "_scib_results.csv")
