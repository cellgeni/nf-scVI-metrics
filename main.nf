#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    ===================
    scVI Hyperparameter Tuning pipeline
    ===================
    This pipeline performs a grid search to optimize hyperparameters for scVI
    Usage: nextflow run main.nf [OPTIONS]
      Required options:
      --input_file: Path to the input file
      Optional options:
      --help: Display this help message
      --umap: Perform UMAP on the scVI embeddings
      --scanvi: Run scANVI label transfer and majority voting aggregation
      --embedding_source: Choose embedding source for scib/umap (scvi or scanvi)
      --seed: Random seed for UMAP and scANVI majority voting
    """.stripIndent()
}

def errorMessage() {
    log.info"""
    ================
    samplefile error
    ================
    You failed to provide the SAMPLEFILE input parameter
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}

def CalculateMemory(adata_size, factor, attempt) {
    if (adata_size < 2.GB) {
        return 2.GB * factor * attempt
    }
    return adata_size * factor * attempt
}

process parse_inputs {
  publishDir 'results', mode: 'copy', pattern: 'input_params.csv'
  input:
    val input_file
  output:
    path 'adata_path', emit: adata_path
    path 'adata_mask_*', emit: adata_mask
    path 'params_*', emit: model_input
    path 'input_params.csv', emit: input_params
  script:
  """
    parse_input.py \
      --input_file '$input_file'
  """
}

process prune_adata {
  memory { CalculateMemory(raw_adata.size(), 4, task.attempt) }
  queue { CalculateMemory(raw_adata.size(), 4, task.attempt) < 680.GB ? "normal" : "hugemem" }
  input:
    path raw_adata
    val input_file
    val adata_mask
  output:
    tuple val(adata_mask), path("pruned_adata_${adata_mask}.h5ad"), emit: adata
    tuple path("pruned_adata_${adata_mask}.h5ad"), path("PCA_params_unintegrated_${adata_mask}.h5ad"), emit: pca
  script:
  """
    prune_adata.py \
      --raw_adata '$raw_adata' \
      --input_file '$input_file' \
      --adata_mask '$adata_mask' \
      --n_pca_unintegrated $params.n_pca_unintegrated
  """
}

process run_scVI {
  publishDir 'results/models', mode: 'copy', pattern: '*.pt'
  memory { CalculateMemory(adata.size(), 4, task.attempt) }
  queue { CalculateMemory(adata.size(), 4, task.attempt) < 680.GB ? "gpu-normal" : "gpu-huge" }
  input:
    tuple val(adata_mask), path(adata), path(model_input)
    val input_file
  output:
    tuple path(adata), path("scvi_${model_input}_${adata_mask}.h5ad"), emit: embedding
    path "history_${model_input}_${adata_mask}", emit: history
    path "model_${model_input}_${adata_mask}.pt", emit: model, optional: true
    tuple val(adata_mask), path(adata), path(model_input), path("model_${model_input}_${adata_mask}.pt"), emit: model_with_inputs, optional: true
  script:
  """
    run_scVI.py \
      --adata '$adata' \
      --input_file '$input_file' \
      --params_file '$model_input' \
      --adata_mask '$adata_mask' \
      --save_model $params.save_model \
      --check_val_every_n_epoch $params.check_val_every_n_epoch
  """
}

process run_scANVI {
  publishDir 'results/scanvi', mode: 'copy', pattern: 'scanvi_*.h5ad'
  memory { CalculateMemory(adata.size(), 4, task.attempt) }
  queue { CalculateMemory(adata.size(), 4, task.attempt) < 680.GB ? "gpu-normal" : "gpu-huge" }
  input:
    tuple val(adata_mask), path(adata), path(model_input), path(scvi_model)
    val input_file
  output:
    path "scanvi_*.h5ad", emit: scanvi_results
    tuple path(adata), path("scanvi_*.h5ad"), emit: scanvi_with_adata
  script:
  """
    run_scANVI.py \
      --adata '$adata' \
      --input_file '$input_file' \
      --scvi_model '$scvi_model' \
      --output_prefix 'scanvi_${scvi_model.baseName}'
  """
}

process scanvi_majority_voting {
  publishDir 'results/scanvi', mode: 'copy'
  input:
    path scanvi_results
  output:
    path 'scanvi_majority.h5ad'
    path 'scanvi_majority.csv'
  script:
  """
    scanvi_majority_voting.py \
      --scanvi_results '$scanvi_results' \
      --output_h5ad scanvi_majority.h5ad \
      --output_csv scanvi_majority.csv \
      --prediction_source $params.scanvi_prediction_source \
      --unlabeled_category $params.scanvi_unlabeled_category \
      --over_clustering_key $params.scanvi_over_clustering_key \
      --resolution $params.scanvi_majority_resolution \
      --seed $params.seed \
      --use_rep $params.scanvi_majority_use_rep
  """
}

process plot_history {
  publishDir 'results', mode: 'copy'
  input:
    path history
  output:
    path 'reconstruction_loss.pdf'
  script:
  """
    plot_history.py \
      --history_files '$history'
  """
}

process run_scib {
  memory { CalculateMemory(adata.size(), 8, task.attempt) }
  queue { CalculateMemory(adata.size(), 8, task.attempt) < 680.GB ? "normal" : "hugemem" }
  input:
    tuple path(adata), path(scVI_embedding)
    val input_file
  output:
    path 'X_*'
  script:
  """
    run_scib.py \
      --adata '$adata' \
      --input_file '$input_file' \
      --embedding_h5ad '$scVI_embedding' \
      --scib_max_obs $params.scib_max_obs \
      --scib_label_key $params.scib_label_key \
      --n_cpu $task.cpus
  """
}

process plot_scib {
  publishDir 'results', mode: 'copy'
  input:
    path results
  output:
    path 'scib_results.svg'
    path 'scib_results_scaled.svg'
  script:
  """
    plot_scib.py \
      --result_files '$results'
  """
}

process run_umap {
  memory { CalculateMemory(adata.size(), 12, task.attempt) }
  queue { CalculateMemory(adata.size(), 12, task.attempt) < 680.GB ? "normal" : "hugemem" }
  input:
    tuple path(adata), path(scVI_embedding)
  output:
    tuple path(adata), path('umap_*')
  script:
  """
    run_umap.py \
      --embedding_h5ad '$scVI_embedding' \
      --seed $params.seed
  """
}

process plot_umap {
  publishDir 'results', mode: 'copy'
  input:
    val input_file
    path umaps
  output:
    path 'umap.pdf'
  script:
  """
    plot_umap.py \
      --input_file '$input_file' \
      --umaps '$umaps'
  """
}

process combine_embedding {
  publishDir 'results', mode: 'copy'
  memory { CalculateMemory(raw_adata.size(), 2, task.attempt) }
  queue { CalculateMemory(raw_adata.size(), 2, task.attempt) < 680.GB ? "normal" : "hugemem" }
  input:
    path raw_adata
    path embeddings
  output:
    path 'adata_embedding.h5ad'
  script:
  """
    comb_embedding.py \
      --raw_adata '$raw_adata' \
      --embeddings '$embeddings'
  """
}

workflow {
  if (params.help) {
    helpMessage()
    exit 0
  }
  else {
    if (params.scanvi && !params.save_model) {
      log.info "ERROR: --scanvi requires --save_model true so scVI models are available for scANVI."
      exit 1
    }
    if (params.embedding_source == 'scanvi' && !params.scanvi) {
      log.info "ERROR: --embedding_source scanvi requires --scanvi true."
      exit 1
    }
    if (params.scib_label_key == 'C_scANVI' && params.embedding_source != 'scanvi') {
      log.info "ERROR: --scib_label_key C_scANVI requires --embedding_source scanvi."
      exit 1
    }
    parse_inputs(params.input_file)
    prune_adata(parse_inputs.out.adata_path.text, params.input_file, parse_inputs.out.adata_mask.text.flatten())

    prune_adata.out.adata
      .combine(parse_inputs.out.model_input.flatten())
      .set {runs}
  
    run_scVI(runs, params.input_file)
    plot_history(run_scVI.out.history.collect())

    run_scVI.out.embedding.concat(prune_adata.out.pca)
      .set {embeddings_scvi}

    if (params.scanvi) {
      run_scANVI(run_scVI.out.model_with_inputs, params.input_file)
      scanvi_majority_voting(run_scANVI.out.scanvi_results.collect())
    }

    if (params.embedding_source == 'scanvi') {
      if (params.scib_label_key == 'C_scANVI') {
        run_scANVI.out.scanvi_with_adata.set {embeddings_scanvi}
      } else {
        run_scANVI.out.scanvi_with_adata.concat(prune_adata.out.pca)
          .set {embeddings_scanvi}
      }
      run_scib(embeddings_scanvi, params.input_file)
      plot_scib(run_scib.out.collect())
      embeddings_only = run_scANVI.out.scanvi_with_adata
      if (params.umap) {
        run_umap(embeddings_scanvi)
      plot_umap(params.input_file, run_umap.out.map { it[1] }.collect())
        embeddings_only = embeddings_only.concat(run_umap.out)
      }
      combine_embedding(parse_inputs.out.adata_path.text, embeddings_only.map { it[1] }.collect())
    } else {
      run_scib(embeddings_scvi, params.input_file)
      plot_scib(run_scib.out.collect())
      embeddings_only = run_scVI.out.embedding
      if (params.umap) {
        run_umap(embeddings_scvi)
        plot_umap(params.input_file, run_umap.out.map { it[1] }.collect())
        embeddings_only = embeddings_only.concat(run_umap.out)
      }
      combine_embedding(parse_inputs.out.adata_path.text, embeddings_only.map { it[1] }.collect())
    }
  }
}
