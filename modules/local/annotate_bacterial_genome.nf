#!/usr/bin/env nextflow

// nextflow.enable.types = true


// use bakta to annotate the genome assembly
process ANNOTATE_BACTERIAL_GENOME {

    container "community.wave.seqera.io/library/bakta:1.12.0--43748ab94e60a85a"

    input:
    tuple val(sampleid), val(taxid), path(assembly_fasta), path(bakta_db)

    output:
    tuple val(sampleid), val(taxid), path("annotation"), emit: all
    tuple val(sampleid), val(taxid), path("annotation/*.fna"), path("annotation/*.gff3"), path("annotation/*.faa"), emit: annotations

    script:
    """
    set -euo pipefail
   
    mkdir -p cache
    mkdir -p config

    export XDG_CACHE_HOME=cache
    export MPLCONFIGDIR=config

    out="annotation"
    mkdir -p "\$out"
   
    bakta --threads ${task.cpus} --force --compliant --output "\$out" --db ${bakta_db} ${assembly_fasta}
    """
}
