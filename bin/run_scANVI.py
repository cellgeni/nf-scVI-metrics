#!/usr/bin/env python3

import argparse
from pathlib import Path
from runpy import run_path
import pandas as pd
import scanpy as sc
import scvi
import torch


def _coerce_labels(series, unlabeled_category: str):
    if unlabeled_category is None:
        return series
    # Replace NaNs with the unlabeled category and cast to string for SCANVI
    if pd.api.types.is_categorical_dtype(series):
        series = series.astype(object)
    series = series.where(~series.isna(), unlabeled_category)
    return series.astype(str)


def main():
    parser = argparse.ArgumentParser(description="Run scANVI label transfer from a trained scVI model.")
    parser.add_argument("--adata", type=str, required=True, help="Path to the pruned AnnData file used for scVI.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the pipeline input params file.")
    parser.add_argument("--scvi_model", type=str, required=True, help="Path to the trained scVI model directory.")
    parser.add_argument("--output_prefix", type=str, default=None, help="Prefix for output files.")
    parser.add_argument("--labels_key", type=str, default=None, help="Override labels key for scANVI.")
    parser.add_argument("--unlabeled_category", type=str, default=None, help="Override unlabeled category.")
    parser.add_argument("--n_samples_per_label", type=int, default=None, help="Override n_samples_per_label for scANVI training.")
    args = parser.parse_args()

    torch.set_float32_matmul_precision("high")

    params = run_path(args.input_file)
    scanvi_input = params.get("scanvi_input", {})
    scanvi_train_input = params.get("scanvi_train_input", {})

    labels_key = args.labels_key or scanvi_input.get("labels_key", "celltype")
    unlabeled_category = args.unlabeled_category or scanvi_input.get("unlabeled_category", "nan")
    n_samples_per_label = args.n_samples_per_label or scanvi_input.get("n_samples_per_label", 100)

    train_kwargs = {
        "max_epochs": scanvi_train_input.get("max_epochs", 30),
        "early_stopping": scanvi_train_input.get("early_stopping", True),
        "early_stopping_patience": scanvi_train_input.get("early_stopping_patience", 20),
        "n_samples_per_label": scanvi_train_input.get("n_samples_per_label", n_samples_per_label),
    }
    # Allow users to pass through any additional train args
    for k, v in scanvi_train_input.items():
        if k not in train_kwargs:
            train_kwargs[k] = v

    adata = sc.read_h5ad(args.adata)
    if labels_key not in adata.obs:
        raise KeyError(f"labels_key '{labels_key}' not found in adata.obs")
    adata.obs[labels_key] = _coerce_labels(adata.obs[labels_key], unlabeled_category)

    scvi_model = scvi.model.SCVI.load(Path(args.scvi_model), adata=adata)

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        unlabeled_category=unlabeled_category,
        labels_key=labels_key,
    )

    scanvi_model.train(**train_kwargs)

    model_base = Path(args.scvi_model).name
    output_prefix = args.output_prefix or f"scanvi_{model_base}"

    scanvi_model.save(f"scanvi_model_{model_base}")

    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)
    adata.obs["C_scANVI"] = scanvi_model.predict(adata)
    adata.uns["df_scANVI"] = scanvi_model.predict(adata, soft=True)

    # Drop expression matrix and layers to keep output small
    if hasattr(adata, "X"):
        del adata.X
    if hasattr(adata, "layers"):
        adata.layers.clear()

    adata.write_h5ad(f"{output_prefix}.h5ad")


if __name__ == "__main__":
    main()
