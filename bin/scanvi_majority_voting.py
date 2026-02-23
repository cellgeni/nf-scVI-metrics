#!/usr/bin/env python3

import argparse
from functools import reduce
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc


def read_soft_predictions(h5ad_path: Path, key: str) -> pd.DataFrame:
    with h5py.File(h5ad_path, "r") as f:
        return ad._io.specs.read_elem(f[f"uns/{key}"])


def read_hard_predictions(h5ad_path: Path, key: str) -> pd.Series:
    with h5py.File(h5ad_path, "r") as f:
        obs = ad._io.specs.read_elem(f["obs"])
    if key not in obs.columns:
        raise KeyError(f"obs key '{key}' not found in {h5ad_path}")
    return obs[key].astype(str)


def build_soft_from_hard(hard_df: pd.DataFrame) -> pd.DataFrame:
    labels = sorted(pd.unique(hard_df.values.ravel()))
    counts = pd.DataFrame(0, index=hard_df.index, columns=labels, dtype=float)
    for col in hard_df.columns:
        counts = counts.add(pd.get_dummies(hard_df[col]), fill_value=0)
    counts = counts.reindex(columns=labels, fill_value=0)
    return counts.div(counts.sum(axis=1).values, axis=0)


def validate_frames(frames):
    if not frames:
        raise ValueError("No scANVI results provided.")
    idx0 = frames[0].index
    col0 = frames[0].columns
    for i, df in enumerate(frames[1:], start=1):
        if not idx0.equals(df.index):
            raise ValueError(f"Index mismatch in scANVI result {i}")
        if not col0.equals(df.columns):
            raise ValueError(f"Column mismatch in scANVI result {i}")


def build_majority_voting(
    adata,
    df_scANVI,
    over_clustering_key,
    resolution,
    seed,
    use_rep,
):
    if over_clustering_key not in adata.obs:
        if use_rep not in adata.obsm:
            raise KeyError(f"use_rep '{use_rep}' not found in adata.obsm")
        sc.pp.neighbors(adata, use_rep=use_rep)
        sc.tl.leiden(
            adata,
            resolution=resolution,
            random_state=seed,
            key_added=over_clustering_key,
            flavor="igraph",
        )

    df_oc_anno = df_scANVI.groupby(adata.obs[over_clustering_key]).mean()
    majority_voting = df_oc_anno.idxmax(axis=1).astype(str).to_dict()
    adata.obs["majority_voting"] = adata.obs[over_clustering_key].map(majority_voting)

    return df_oc_anno




def main():
    parser = argparse.ArgumentParser(description="Aggregate scANVI soft predictions and compute majority voting labels.")
    parser.add_argument("--scanvi_results", required=True, help="Paths to scANVI result h5ad files.")
    parser.add_argument("--output_h5ad", required=True, help="Path to output h5ad file.")
    parser.add_argument("--output_csv", required=True, help="Path to output CSV file.")
    parser.add_argument("--soft_key", default="df_scANVI", help="Key in adata.uns for soft predictions.")
    parser.add_argument("--hard_key", default="C_scANVI", help="Key in adata.obs for hard predictions.")
    parser.add_argument("--prediction_source", choices=["soft", "hard"], default="soft", help="Use soft (df_scANVI) or hard (C_scANVI) predictions.")
    parser.add_argument("--unlabeled_category", default="nan", help="Label used for unlabeled cells.")
    parser.add_argument("--over_clustering_key", default="over_clustering", help="obs key for over-clustering.")
    parser.add_argument("--resolution", type=float, default=30.0, help="Leiden resolution for over-clustering.")
    parser.add_argument("--seed", type=int, default=123, help="Random seed for Leiden.")
    parser.add_argument("--use_rep", default="X_scANVI", help="Embedding key in adata.obsm.")
    args = parser.parse_args()

    scanvi_results = [Path(p) for p in args.scanvi_results.split()]
    adata = sc.read_h5ad(scanvi_results[0])

    if args.prediction_source == "soft":
        soft_frames = [read_soft_predictions(p, args.soft_key) for p in scanvi_results]
        validate_frames(soft_frames)

        stacked = np.stack([df.to_numpy() for df in soft_frames])
        mean_soft = stacked.mean(axis=0)
        df_scANVI = pd.DataFrame(mean_soft, index=soft_frames[0].index, columns=soft_frames[0].columns)
    else:
        hard_series = [read_hard_predictions(p, args.hard_key) for p in scanvi_results]
        hard_df = pd.concat(hard_series, axis=1)
        hard_df.columns = [p.stem for p in scanvi_results]
        df_scANVI = build_soft_from_hard(hard_df)

    adata.uns[args.soft_key] = df_scANVI

    majority = df_scANVI.idxmax(axis=1).astype(str)
    adata.obs["majority"] = majority.reindex(adata.obs_names).values

    df_oc_anno = build_majority_voting(
        adata,
        df_scANVI,
        args.over_clustering_key,
        args.resolution,
        args.seed,
        args.use_rep,
    )


    df_oc_anno = df_oc_anno.rename(columns={i: f"mv_{i}" for i in df_oc_anno.columns})
    adata.uns["df_oc_anno"] = df_oc_anno.copy()

    comb_df = pd.concat(
        [
            adata.obs[["majority", args.over_clustering_key, "majority_voting"]],
            df_scANVI,
        ],
        axis=1,
    )

    comb_df = comb_df.merge(df_oc_anno, left_on=args.over_clustering_key, right_index=True, how="left")

    adata.write_h5ad(args.output_h5ad)
    comb_df.to_csv(args.output_csv)


if __name__ == "__main__":
    main()
