#!/usr/bin/env python3

import argparse
import json
import pandas as pd
from runpy import run_path

parser = argparse.ArgumentParser(description='Parse input parameters.')
parser.add_argument('--input_file', type=str, help='Path to the parameters file.')
args = parser.parse_args()

params = run_path(args.input_file)

def create_turning_grid(**kwargs):
    out = [{}]
    for in_key, in_val in kwargs.items():
        temp = out
        out = []
        for elem in temp:
            for in_x in in_val:
                temp_elem = elem.copy()
                temp_elem.update({in_key: in_x})
                out.append(temp_elem)
    return out

with open("adata_path", "w") as f:
    f.write(params['anndata_input'])

turning_grid = create_turning_grid(**params['model_input'])

out_dict = dict()
for i, j in enumerate(turning_grid):
    with open(f"params_{i}", "w") as f:
        json.dump(j, f)
    out_dict[f"params_{i}"] = j

pd.DataFrame.from_dict(out_dict).transpose().to_csv('input_params.csv')
