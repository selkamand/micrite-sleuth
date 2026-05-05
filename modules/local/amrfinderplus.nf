#!/usr/bin/env nextflow

process AMRFINDER {
    tag "${sampleid}:${taxid}"
    container "community.wave.seqera.io/library/bakta:1.12.0--43748ab94e60a85a"
    cpus 1

    input:
    tuple val(sampleid), val(taxid), path(fna), path(faa), path(gff3)
    path database

    output:
    tuple val(sampleid), val(taxid), path("amrfinder.tsv")

    script:
    """
    set -euo pipefail

    amrfinder --plus --database ${database} -n ${fna} -p ${faa} --gff ${gff3} --annotation_format bakta > amrfinder.tsv
    """
}
