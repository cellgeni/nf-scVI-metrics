#!/usr/bin/env python3

import sys
import pandas as pd
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

result_list = []
for i in sorted(sys.argv[1].split()):
    result_list.append(pd.read_csv(i, index_col = 0).iloc[:,:-1])
result_list.append(pd.read_csv(i, index_col = 0).iloc[:,-1:])

bm = Benchmarker(None, None, None, [None], 
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection())

bm._results = pd.concat(result_list, axis = 1)
bm._embedding_obsm_keys = sys.argv[1].split()
bm._benchmarked = True

bm.plot_results_table(show = False).figure.savefig('scib_results_scaled.svg')
bm.plot_results_table(min_max_scale = False, show = False).figure.savefig('scib_results.svg')