#!/usr/bin/env nextflow

// Use krakentools to extract all reads classified as the given taxid (including children) 
// krakentools is single threaded, so don't bother wasting extra cores on this step.
process EXTRACT_READS_BY_TAXID {
    cpus 1

    input:
    tuple val(sampleid), val(taxid), path(kraken), path(kreport), path(r1), path(r2)

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.taxid_${taxid}.R1.fq"), path("${sampleid}.taxid_${taxid}.R2.fq")

    script:
    """
    extract_kraken_reads.py \
    --include-children -k ${kraken} \
    --fastq-output \
    -r ${kreport} \
    -t ${taxid} \
    -s ${r1} -s2 ${r2} \
    -o ${sampleid}.taxid_${taxid}.R1.fq -o2 ${sampleid}.taxid_${taxid}.R2.fq 
  """
}
