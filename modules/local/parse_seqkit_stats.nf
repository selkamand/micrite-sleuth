#!/usr/bin/env nextflow

// Get total  to get gene completeness metrics 
process COUNT_TOTAL_BASES {
    input:
    tuple val(sampleid), val(taxid), path(seqkit_stats)

    output:
    tuple val(sampleid), val(taxid), path("total_length.txt")

    script:
    """
    set -euo pipefail

    awk 'NR>1 { val += \$5 } END { print val }' ${seqkit_stats} > total_length.txt
    """
}
