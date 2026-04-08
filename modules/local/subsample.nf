#!/usr/bin/env nextflow

// Use fq to substample nreads then convert to fasta to run blast
process SUBSAMPLE_FQ {
    ext subsample_blast_seed: 111

    input:
    tuple val(sampleid), val(taxid), path(r1), val(nreads)

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.taxid_${taxid}.subsampled.${nreads}reads.fasta")

    script:
    """
    set -euo pipefail

    # Subsample
    fq subsample \
        --record-count ${nreads} \
        --seed ${task.ext.subsample_blast_seed} \
        --r1-dst ${sampleid}.taxid_${taxid}.subsampled.${nreads}reads.fq \
        ${r1}

    # Convert to fasta
    seqkit fq2fa ${sampleid}.taxid_${taxid}.subsampled.${nreads}reads.fq > ${sampleid}.taxid_${taxid}.subsampled.${nreads}reads.fasta
    """
}
