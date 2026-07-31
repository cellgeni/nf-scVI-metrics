#!/usr/bin/env python3

import argparse
import itertools
import json
from pathlib import Path

import pandas as pd

try:
    import yaml
except ImportError as exc:
    raise ImportError(
        "PyYAML is required to parse YAML input files. Please install pyyaml."
    ) from exc


def load_yaml_input(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError("Input YAML must contain a top-level mapping/object.")
    return data


def validate_params(params):
    required_keys = ["preprocess", "scvi", "scib"]
    for key in required_keys:
        if key not in params:
            raise KeyError(f"Missing required key in input YAML: '{key}'")

    preprocess = params["preprocess"]
    scvi_config = params["scvi"]
    scib_config = params["scib"]

    if not isinstance(preprocess, dict):
        raise ValueError("'preprocess' must be a mapping/object.")
    if not isinstance(scvi_config, dict):
        raise ValueError("'scvi' must be a mapping/object.")
    if not isinstance(scib_config, dict):
        raise ValueError("'scib' must be a mapping/object.")

    if "anndata_input" not in preprocess:
        raise KeyError("Missing required key in input YAML: 'preprocess.anndata_input'")
    if not isinstance(preprocess["anndata_input"], str) or not preprocess["anndata_input"].strip():
        raise ValueError("'anndata_input' must be a non-empty string.")

    if "model_input" not in scvi_config:
        raise KeyError("Missing required key in input YAML: 'scvi.model_input'")
    model_input = scvi_config["model_input"]
    if not isinstance(model_input, dict) or not model_input:
        raise ValueError("'scvi.model_input' must be a non-empty mapping of parameter lists.")

    for param_name, values in model_input.items():
        if not isinstance(values, list) or not values:
            raise ValueError(
                f"'scvi.model_input.{param_name}' must be a non-empty list."
            )

    anndata_mask = preprocess.get("anndata_mask", [])
    if not isinstance(anndata_mask, list):
        raise ValueError("'preprocess.anndata_mask' must be a list.")

    if "setup_anndata" not in scvi_config or not isinstance(scvi_config["setup_anndata"], dict):
        raise ValueError("'scvi.setup_anndata' must be provided as a mapping/object.")
    if "train" not in scvi_config or not isinstance(scvi_config["train"], dict):
        raise ValueError("'scvi.train' must be provided as a mapping/object.")
    if "metrics" not in scib_config or not isinstance(scib_config["metrics"], dict):
        raise ValueError("'scib.metrics' must be provided as a mapping/object.")


def create_tuning_grid(model_input):
    keys = list(model_input.keys())
    value_product = itertools.product(*(model_input[key] for key in keys))
    return [dict(zip(keys, values)) for values in value_product]


def main():
    parser = argparse.ArgumentParser(description="Parse YAML input parameters.")
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to the YAML parameters file.",
    )
    args = parser.parse_args()

    params = load_yaml_input(args.input_file)
    validate_params(params)

    with open("adata_path", "w", encoding="utf-8") as handle:
        handle.write(params["preprocess"]["anndata_input"])

    tuning_grid = create_tuning_grid(params["scvi"]["model_input"])

    out_dict = {}
    for index, item in enumerate(tuning_grid):
        with open(f"params_{index}", "w", encoding="utf-8") as handle:
            json.dump(item, handle)
        out_dict[f"params_{index}"] = item

    pd.DataFrame.from_dict(out_dict).transpose().to_csv("input_params.csv")

    masks = params["preprocess"].get("anndata_mask", [])
    if not masks:
        masks = [""]
    for index, mask in enumerate(masks):
        with open(f"adata_mask_{index}", "w", encoding="utf-8") as handle:
            handle.write(str(mask))


if __name__ == "__main__":
    main()