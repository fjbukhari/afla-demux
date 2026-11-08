#!/usr/bin/env nextflow
// AFLA demultiplexer: sorts Illumina or Nanopore reads into one FASTQ set per sample from a sample sheet.
// Becomes main.nf in the shared package. Education and research use only, not for clinical use.

process DEMUX_ILLUMINA {
    container "quay.io/biocontainers/python:3.12"
    cpus 2
    publishDir "${params.out_dir}", mode: "copy"
    input:
    path input
    path sheet
    output:
    path "demux/**"
    script:
    """
    afla_demux_illumina.py --sample-sheet ${sheet} --input ${input} --out demux \\
      --mismatches ${params.mismatches} --i5 ${params.i5_direction}
    """
}

process DEMUX_ONT {
    label "big"
    container "quay.io/biocontainers/python-edlib:1.3.9.post1--py312h248b54c_3"
    publishDir "${params.out_dir}", mode: "copy"
    input:
    path input
    path sheet
    output:
    path "demux/**"
    script:
    def both = params.require_both_ends.toString() == 'true' ? '--require-both' : ''
    def maxlen = params.max_length ? "--max-length ${params.max_length}" : ''
    """
    afla_demux.py --sample-sheet ${sheet} --input ${input} --out demux --threads ${task.cpus} \\
      --max-error ${params.max_error} --min-margin ${params.min_margin} --window ${params.search_window} \\
      --min-length ${params.min_length} ${maxlen} ${both}
    """
}

workflow {
    if (!(params.platform in ['illumina', 'ont'])) {
        error "platform must be 'illumina' or 'ont'."
    }
    if (!params.input || !params.sample_sheet) {
        error "Give the folder of reads (input) and the sample sheet."
    }
    def input = file(params.input, checkIfExists: true)
    def sheet = file(params.sample_sheet, checkIfExists: true)
    if (params.platform == 'illumina') {
        DEMUX_ILLUMINA(input, sheet)
    } else {
        DEMUX_ONT(input, sheet)
    }
}
