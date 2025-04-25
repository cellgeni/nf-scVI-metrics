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
    result = [[]]
    for in_val in kwargs.values():
        result = [x+[y] for x in result for y in in_val]
    return [dict(zip(kwargs.keys(), i)) for i in result]

with open("adata_path", "w") as f:
    f.write(params['anndata_input'])

turning_grid = create_turning_grid(**params['model_input'])

out_dict = dict()
for i, j in enumerate(turning_grid):
    with open(f"params_{i}", "w") as f:
        json.dump(j, f)
    out_dict[f"params_{i}"] = j

pd.DataFrame.from_dict(out_dict).transpose().to_csv('input_params.csv')

temp = params['anndata_mask']
if temp == []:
    temp = ['']
for i, j in enumerate(temp):
    with open(f"adata_mask_{i}", "w") as f:
        f.write(j)