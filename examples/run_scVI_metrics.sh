module load cellgen/nextflow
module load cellgen/singularity

nextflow run main.nf \
    --input_file '/lustre/scratch127/cellgen/cellgeni/nf-scVI-metrics/examples/inputs.py' \
    --umap \
    -resume

