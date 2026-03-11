#!/usr/bin/env nexftlow

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
