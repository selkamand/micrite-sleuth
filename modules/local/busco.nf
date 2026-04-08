#!/usr/bin/env nextflow

// Use busco to get  completeness metrics 
process BUSCO_COMPLETENESS {

    container "community.wave.seqera.io/library/busco:6.0.0--a9a1426105f81165"

    input:
    tuple val(sampleid), val(taxid), path(assembly_fasta)
    // A valid busco lineage (fed to --lineage argument). See busco --list-datasets for all options
    val lineage
    // datasets is a directory that contains the busco lineage datasets. 
    path datasets

    output:
    tuple val(sampleid), val(taxid), path("busco.${lineage}.txt")

    script:
    """
    set -euo pipefail

    busco \
    --offline \
    --download_path ${datasets} \
    --opt-out-run-stats \
    -m genome \
    -l ${lineage} \
    --in ${assembly_fasta} \
    -o busco.${lineage}.txt
    """
}
