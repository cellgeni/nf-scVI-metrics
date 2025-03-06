import sys
import json
import pandas as pd

temp = dict()
for i in sorted(sys.argv[0].split()):
    with open(i) as f:
        model_input = json.load(f)
    temp[i] = model_input
pd.DataFrame.from_dict(temp).transpose().to_csv('input_params.csv')
