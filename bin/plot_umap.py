#!/usr/bin/env python3

import argparse
import anndata as ad
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from runpy import run_path

parser = argparse.ArgumentParser(description='Plot UMAP results.')
parser.add_argument('--input_file', type=str, help='Path to the parameters file.')
parser.add_argument('--umaps', type=str, help='Space-separated list of UMAP h5ad files.')
args = parser.parse_args()

params = run_path(args.input_file)

pp = PdfPages("umap.pdf")
sc.set_figure_params(figsize=(8, 8))
for i in sorted(args.umaps.split()):
    adata = ad.read_h5ad(i, backed='r')
    label_key = params['scib_input']['label_key']
    if label_key not in adata.obs:
        if 'C_scANVI' in adata.obs:
            label_key = 'C_scANVI'
        else:
            raise KeyError(f"label_key '{label_key}' not found in {i} obs.")
    plt.clf()
    fig = sc.pl.umap(adata, color=[params['scib_input']['batch_key'], label_key], 
                     return_fig=True, show=False)
    fig.suptitle(i.split('.')[0])
    pp.savefig(fig, bbox_inches='tight')
pp.close()
