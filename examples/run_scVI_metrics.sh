module load cellgen/nextflow
module load cellgen/singularity

nextflow run ~/nf-scVI-metrics/main.nf \
    --input_file '~/nf-scVI-metrics/examples/inputs.py' \
    --umap \
    --save_model \
    --scanvi \
    --embedding_source scanvi \
    --scib_label_key C_scANVI \
    -resume
    # -profile local \
