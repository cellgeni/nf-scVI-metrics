import sys
import h5py
import anndata as ad
import scanpy as sc

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from runpy import run_path


adata = ad.read_h5ad(sys.argv[1], backed = True)
args = run_path(sys.argv[2])

pp = PdfPages("umap.pdf")
sc.set_figure_params(figsize=(8, 8))
for i in sorted(sys.argv[3].split()):
    with h5py.File(i) as file:
        umap = ad._io.specs.read_elem(file['obsm/X_umap'])
    adata.obsm['X_umap'] = umap
    plt.clf()
    fig = sc.pl.umap(adata, color = [args['scib_input']['batch_key'], args['scib_input']['label_key']], 
                     return_fig = True, show = False)
    fig.suptitle(i)
    pp.savefig(fig, bbox_inches='tight')
pp.close()
