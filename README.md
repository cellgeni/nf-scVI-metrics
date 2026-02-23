# scVI Hyperparameter Metrics Pipeline  

This pipeline performs a grid search to optimize hyperparameters for `scVI`. The grid is constructed by combining elements from predefined lists of model input parameters. After running `scVI` with each set of parameters, integration metrics are calculated using `scib-metrics`.  

## Repository Contents  

- **`main.nf`** – The Nextflow pipeline that orchestrates the entire workflow.  
- **`nextflow.config`** – Configuration file that:  
  - Submits processes to LSF on Sanger's HPC.  
  - Ensures the correct environment is set via a Singularity container (absolute path).  
  - Defines global default parameters.  
- **`examples/inputs.py`** – Example `.py` file specifying the path to `anndata` and configurations for `scVI` and `scib-metrics`.  
- **`examples/run_scVI_metrics.sh`** – Example script to execute the pipeline. 
- **`Dockerfile`** – Defines a Docker image with `scVI` and `scib-metrics` version `0.2.2`.  
- **`bin/run_scANVI.py`** – Runs scANVI label transfer using a trained scVI model.  
- **`bin/scanvi_majority_voting.py`** – Aggregates scANVI soft predictions and computes majority voting labels.  
- **Embedding outputs** – scVI/scANVI embeddings are exported as `.h5ad` with `X_scVI`/`X_scANVI` in `.obsm`.  

## Pipeline Arguments  

- **`--input_file`** – Path to a `.py` file specifying:  
  - The path to `anndata`.  
  - Configuration details for `scVI` and `scib-metrics`. 
  - Optional `scvi_input`, `scvi_train_input`, `scanvi_input`, and `scanvi_train_input` blocks.

### Optional parameters:
* `--help` — Display this help message
* `--umap` — Set this flag to calculate umap
* `--scanvi` — Run scANVI label transfer and majority voting aggregation (requires `--save_model true`)
* `--scanvi_prediction_source` — Use `soft` (df_scANVI) or `hard` (C_scANVI) predictions for majority voting (default: `soft`)
* `--embedding_source` — Choose embedding source for scIB/UMAP (`scvi` or `scanvi`, default: `scvi`)
* `--seed` — Random seed for UMAP and scANVI majority voting (default: `123`)
* `--scib_label_key` — Override `scib_input.label_key` (e.g., `C_scANVI`)
