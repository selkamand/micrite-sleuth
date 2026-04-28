#!/usr/bin/env nextflow

// use barrnap to detect RNA features (inc. 16S/23S rRNA) from de novo assembly 
process BARRNAP {

    container "community.wave.seqera.io/library/barrnap:0.9--3a151f9727225d80"

    input:
    tuple val(sampleid), val(taxid), path(assembly_fasta)

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.${taxid}.barrnap.fasta")

    script:
    """
    set -euo pipefail
     
    barrnap --outseq ${sampleid}.${taxid}.barrnap.fasta ${assembly_fasta}
    """
}


// Split barrnap fasta into 3 fastqs with a single entry (16S SSU, 23S LSU & 5S LSU) using seqkit.
// This lets us use databases specific to each for downstream classification (and multiple sequence alignment)
process SPLIT_RRNA_FASTA {
    cpus 1

    input:
    tuple val(sampleid), val(taxid), path(barrnap_fasta)

    output:
    tuple val(sampleid), val(taxid), path("${barrnap_fasta.baseName}.16S.fasta"), emit: SSU_16S
    tuple val(sampleid), val(taxid), path("${barrnap_fasta.baseName}.23S.fasta"), emit: LSU_23S
    tuple val(sampleid), val(taxid), path("${barrnap_fasta.baseName}.5S.fasta"), emit: LSU_5S

    script:
    """
    set -euo pipefail

    # Grab the first barrnap fasta entry that starts with 16S
    # We do this because sometimes barrnap pulls both a + and -ve strand hit.
    seqkit grep --use-regexp -p ^16S ${barrnap_fasta} | seqkit head -n 1 > ${barrnap_fasta.baseName}.16S.fasta
    seqkit grep --use-regexp -p ^23S ${barrnap_fasta} | seqkit head -n 1 > ${barrnap_fasta.baseName}.23S.fasta
    seqkit grep --use-regexp -p ^5S ${barrnap_fasta} | seqkit head -n 1 > ${barrnap_fasta.baseName}.5S.fasta
    """
}
