#!/usr/bin/env nextflow

// Use fq to substample nreads then convert to fasta to run blast
process SUBSAMPLE_BY_PROPORTION {
    ext seed: 999

    input:
    tuple val(sampleid), val(taxid), path(r1), path(r2), val(subsample_prop)

    output:
    tuple val(sampleid), val(taxid), path("subsampled.R1.fq"), path("subsampled.R2.fq")

    script:
    """
    set -euo pipefail

    # Subsample
    fq subsample --probability ${subsample_prop} --seed ${task.ext.seed} \
    --r1-dst subsampled.R1.fq \
    --r2-dst subsampled.R2.fq \
    ${r1} ${r2}

    """
}
