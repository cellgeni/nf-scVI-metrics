# path to anndata object
anndata_input = "/lustre/scratch127/cellgen/cellgeni/GBM_scVI_metrics.h5ad"

# parameters for scvi.model.SCVI
# please put all parameters as a list
model_input = {
    'n_hidden': [64, 128], 
    'n_latent': [5, 10, 20], 
    'n_layers': [1, 2], 
    'dropout_rate': [0.1, 0.2]
}

# parameters for SCVI.setup_anndata
# adata does no need to be specified
param_input = {
#    'layer': "counts",
    'batch_key': "donor_id",
    'continuous_covariate_keys': ["log1p_n_genes_by_counts", "pct_counts_mt"],
    'categorical_covariate_keys': ["sample", "data_source", "dataset", "phase"]
}

# parameters for scib_metrics.benchmark.Benchmarker
scib_input = {
    'batch_key': "donor_id",
    'label_key': "annotation_coarse"
}

# parameters for SCVI.train
train_input = {
    'train_size': 0.95, 
    'max_epochs': 200, 
    'batch_size': 512
}
