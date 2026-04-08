#!/usr/bin/env nextflow

// Classify ssRNA (barrnap output) using the SILVA database ACT approach 
// Use SINA to align our extracted rRNA sequences to a SILVA SSU (16S) database 
process SINA_ALIGN {

    container "community.wave.seqera.io/library/sina:1.7.2--322a08ea99ba083b"

    input:
    tuple val(sampleid), val(taxid), path(fasta_16s)
    path silva_arb_database

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.${taxid}.msa.fasta")

    script:
    """
    set -euo pipefail
    sina -i ${fasta_16s} -r ${silva_arb_database} -o ${sampleid}.${taxid}.msa.fasta 
    # Consider adding the --search option to classify the result
    """
}
