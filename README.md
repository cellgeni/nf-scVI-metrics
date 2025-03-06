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

## Pipeline Arguments  

- **`--input_file`** – Path to a `.py` file specifying:  
  - The path to `anndata`.  
  - Configuration details for `scVI` and `scib-metrics`. 

### Optional parameters:
* `--help` — Display this help message
* `--umap` — Set this flag to calculate umap
