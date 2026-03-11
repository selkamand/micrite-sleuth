#!/usr/bin/env nexftlow

// Use quast to get whole genome assembly stats
process QUAST_WHOLE_GENOME_ASSEMBLY {

    container "selkamandcci/quast:5.2.0"

    input:
    tuple val(sampleid), val(taxid), path(assembly_fasta), val(min_contig), val(est_ref_size)

    output:
    tuple val(sampleid), val(taxid), path("quast")

    script:
    """
    set -euo pipefail

    quast.py -o quast \
        --min-contig ${min_contig} \
        --est-ref-size ${est_ref_size} \
        ${assembly_fasta}
    """
}
