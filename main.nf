#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    ===================
    scVI Hyperparameter Tuning pipeline
    ===================
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

process parse_inputs {
  publishDir 'results', mode: 'copy', pattern: 'input_params.csv'
  input:
  	val input_file
  output:
	  path 'adata_path', emit: adata_path
    path 'params_*', emit: model_input
    path 'input_params.csv', emit: input_params
  script:
  """
	  parse_input.py '$input_file'
  """
}

process prune_adata {
  memory {raw_adata.size() < 2.GB ? 6.GB * task.attempt : raw_adata.size() * 3 * task.attempt}
  input:
    path raw_adata
    val input_file
  output:
    path 'input_adata', emit: input_adata
    path 'PCA_params_unintegrated.npy', emit: pca
  script:
  """
	  prune_adata.py '$raw_adata' '$input_file'
  """
}

process run_scVI {
  memory {input_adata.size() < 2.GB ? 8.GB * task.attempt : input_adata.size() * 4 * task.attempt}
  input:
    path input_adata
    val input_file
    path model_input
  output:
    path "scvi_${model_input}.npy", emit: embedding
    path "history_${model_input}", emit: history
  script:
  """
	  run_scVI.py '$input_adata' '$input_file' '$model_input'
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
	  plot_history.py '$history'
  """
}

process run_scib {
  memory {input_adata.size() < 2.GB ? 16.GB * task.attempt : input_adata.size() * 8 * task.attempt}
  input:
    path input_adata
    val input_file
    path embedding
  output:
	  path 'param_*'
  script:
  """
    run_scib.py '$input_adata' '$input_file' '$embedding'
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
    plot_scib.py '$results'
  """
}

process run_umap {
  memory {input_adata.size() < 2.GB ? 16.GB * task.attempt : input_adata.size() * 8 * task.attempt}
  input:
    path input_adata
    path scVI_embedding
  output:
    path 'umap_*' 
  script:
  """
	  run_umap.py '$scVI_embedding'
  """
}

process plot_umap {
  publishDir 'results', mode: 'copy'
  input:
	  path input_adata
    val input_file
    path umaps
  output:
    path 'umap.pdf'
  script:
  """
    plot_umap.py '$input_adata' '$input_file' '$umaps'
  """
}

process combine_embedding {
  publishDir 'results', mode: 'copy'
  memory {raw_adata.size() < 2.GB ? 4.GB * task.attempt : raw_adata.size() * 2 * task.attempt}
  input:
    path raw_adata
    path embeddings
  output:
	  path 'adata_embedding.h5ad'
  script:
  """
    comb_embedding.py '$raw_adata' '$embeddings'
  """
}

workflow {
  if (params.help) {
    helpMessage()
    exit 0
  }
  else {
	  parse_inputs(params.input_file)
	  prune_adata(parse_inputs.out.adata_path.text, params.input_file)
	  run_scVI(prune_adata.out.input_adata, params.input_file, parse_inputs.out.model_input.flatten())
	  plot_history(run_scVI.out.history.collect())
	  run_scib(prune_adata.out.input_adata, params.input_file, run_scVI.out.embedding.concat(prune_adata.out.pca))
	  plot_scib(run_scib.out.collect())
    embeddings = run_scVI.out.embedding
	  if (params.umap) {
	    run_umap(prune_adata.out.input_adata, run_scVI.out.embedding.concat(prune_adata.out.pca))
      plot_umap(prune_adata.out.input_adata, params.input_file, run_umap.out.collect())
      embeddings.concat(run_umap.out)
      combine_embedding(parse_inputs.out.adata_path.text, embeddings.collect())
	  } else {
      combine_embedding(parse_inputs.out.adata_path.text, embeddings.collect())
    }
  }
}


