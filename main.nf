#!/usr/bin/env nextflow

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
	path 'model_input_*', emit: model_input
	path 'input_params.csv', emit: input_params
  script:
  """
	python ${projectDir}/bin/parse_input.py '$input_file'
  """
}


process get_adata {
  input:
        path adata_path
  output:
        path 'raw_adata'
  script:
  """
	cp \$(cat ${'adata_path'}) 'raw_adata'
  """
}

process prune_adata {
  memory {raw_adata.size() < 2.GB ? 4.GB * task.attempt : raw_adata.size() * 2 * task.attempt}
  input:
	path raw_adata
	val input_file
  output:
	path 'input_adata'
  script:
  """
	python ${projectDir}/bin/prune_adata.py '$raw_adata' '$input_file'
  """
}

process run_scVI {
  memory {input_adata.size() < 2.GB ? 8.GB * task.attempt : input_adata.size() * 4 * task.attempt}
  input:
	path input_adata
	val input_file
	path model_input
  output:
        path "adata_${model_input}", emit: adata
	path "history_${model_input}", emit: history
  script:
  """
	python ${projectDir}/bin/run_scVI.py '$input_adata' '$input_file' '$model_input'
  """
}

process run_umap {
  memory {input_adata.size() < 2.GB ? 16.GB * task.attempt : input_adata.size() * 8 * task.attempt}
  input:
	path input_adata
	path adata_scVI
  output:
	path 'umap_*' 
  script:
  """
	python ${projectDir}/bin/run_umap.py '$adata_scVI'
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
	python ${projectDir}/bin/comb_history.py '$history'
  """
}

process run_scib {
  memory {input_adata.size() < 2.GB ? 32.GB * task.attempt : input_adata.size() * 32 * task.attempt}
  publishDir 'results', mode: 'copy'
  input:
        path input_adata
        val input_file
        path adata
  output:
	path 'adata_embedding.h5ad'
        path 'scib_results.svg'
	path 'scib_results_scaled.svg'
  script:
  """
        python ${projectDir}/bin/run_scib.py '$input_adata' '$input_file' '$adata'
  """
}

process plot_umap {
  publishDir 'results', mode: 'copy'
  input:
	path input_adata
        val input_file
        path adata
  output:
        path 'umap.pdf'
  script:
  """
        python ${projectDir}/bin/plot_umap.py '$input_adata' '$input_file' '$adata'
  """
}

workflow {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
	parse_inputs(params.input_file)
	get_adata(parse_inputs.out.adata_path)
	prune_adata(get_adata.out, params.input_file)
	run_scVI(prune_adata.out, params.input_file, parse_inputs.out.model_input.flatten())
	plot_history(run_scVI.out.history.collect())
	run_scib(prune_adata.out, params.input_file, run_scVI.out.adata.collect())
	if (params.umap) {
	  run_umap(prune_adata.out, run_scVI.out.adata)
          plot_umap(prune_adata.out, params.input_file, run_umap.out.collect())
	}
  }
}


