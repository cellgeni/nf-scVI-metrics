import sys
import h5py
import anndata as ad
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from runpy import run_path

adata = ad.read_h5ad(sys.argv[1])
args = run_path(sys.argv[2])
temp = sorted(sys.argv[3].split())

if 'layer' in args['param_input']:
        adata.X = adata.layers[args['param_input']['layer']]

list_obsm = []
from runpy import run_path
for i in temp:
    with h5py.File(i) as file:
        embedding = ad._io.specs.read_elem(file['obsm/scVI'])
    temp_obsm = 'param_' + i.split('_')[3]
    adata.obsm[temp_obsm] = embedding
    list_obsm.append(temp_obsm)

bm = Benchmarker(
    adata,
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    embedding_obsm_keys=list_obsm,
    n_jobs=8,
    **args['scib_input']
)


adata.write_h5ad("adata_embedding.h5ad")

bm.benchmark()

bm.plot_results_table(show = False).figure.savefig('scib_results_scaled.svg')
bm.plot_results_table(min_max_scale = False, show = False).figure.savefig('scib_results.svg')
