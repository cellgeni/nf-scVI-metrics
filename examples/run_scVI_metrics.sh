module load cellgen/nextflow
module load cellgen/singularity

nextflow run main.nf --input_file '/lustre/scratch127/cellgen/cellgeni/yc6/Hyperparameter_scVI/nf-scHPT/examples/inputs.py' -resume

