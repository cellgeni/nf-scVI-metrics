#!/usr/bin/env python3

import sys
import argparse
import h5py
import anndata as ad
import scanpy as sc
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from runpy import run_path

parser = argparse.ArgumentParser(description='Plot UMAP results.')
parser.add_argument('--adata', type=str, help='Path to the AnnData file.')
parser.add_argument('--input_file', type=str, help='Path to the parameters file.')
parser.add_argument('--umaps', type=str, help='Space-separated list of UMAP numpy files.')
args = parser.parse_args()

adata = ad.read_h5ad(args.adata, backed=True)
params = run_path(args.input_file)

pp = PdfPages("umap.pdf")
sc.set_figure_params(figsize=(8, 8))
for i in sorted(args.umaps.split()):
    umap = np.load(i)
    adata.obsm['X_umap'] = umap
    plt.clf()
    fig = sc.pl.umap(adata, color=[params['scib_input']['batch_key'], params['scib_input']['label_key']], 
                     return_fig=True, show=False)
    fig.suptitle(i.split('.')[0])
    pp.savefig(fig, bbox_inches='tight')
pp.close()
