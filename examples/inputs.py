# path to anndata object
anndata_input = "/nfs/cellgeni/yc6/test_adata.h5ad"

# List of sets of genes to use for scVI. Use strings referring to arrays in var.
# '' means all genes, use all genes by default([]) 
anndata_mask = [
    'hvg_ds_3k', 
    # 'hvg_ds_5k',
    'hvg_sid_3k',
    # 'hvg_sid_5k',
    ]

# parameters for scvi.model.SCVI
# please put all parameters as a list
model_input = {
    'n_hidden': [128], 
    'n_latent': [10, 20], 
    'n_layers': [1], 
    'dispersion': ['gene'],
    'gene_likelihood': ['nb'],
    'use_observed_lib_size': [False],
    # 'dropout_rate': [0.1, 0.2]
}

# parameters for SCVI.setup_anndata
# adata does no need to be specified
scvi_input = {
    # 'layer': "counts",
    'batch_key': "sample_id",
    'continuous_covariate_keys': ["total_counts", "n_genes_by_counts", "pct_counts_MT"],
    'categorical_covariate_keys': ["donor_id", "tissue_type", "kit_10X", "phase"]
}

# parameters for scib_metrics.benchmark.Benchmarker
scib_input = {
    'batch_key': "sample_id",
    'label_key': "C_scANVI",
}

# parameters for SCVI.train
scvi_train_input = {
    'train_size': 0.95, 
    'max_epochs': 50, 
    'batch_size': 512,
    'early_stopping': True,
    'early_stopping_patience': 10,
}

# parameters for scANVI label transfer
scanvi_input = {
    'labels_key': "celltype",
    'unlabeled_category': "nan",
    'n_samples_per_label': 10
}

# parameters for scANVI.train
scanvi_train_input = {
    'max_epochs': 30,
    'early_stopping': True,
    'early_stopping_patience': 20,
}
