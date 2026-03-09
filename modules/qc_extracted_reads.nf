#!/usr/bin/env nexftlow

// Run FASTQC and seqkit stats on a pair of reads belonging to a particular taxid. 
// Assumes reads ahave been extracted by EXTRACT_READS_BY_TAXID

process QC_EXTRACTED_READS {
    tag "${sampleid}.${taxid}"

    input:
    tuple val(sampleid), val(taxid), path(fq1), path(fq2)

    output:
    tuple val(sampleid), val(taxid), path("read_stats")

    script:
    """
    outdir="read_stats"
    mkdir -p "\${outdir}"
    fastqc --nogroup -o "\${outdir}" ${fq1} ${fq2}
    seqkit stats ${fq1} ${fq2} > "\${outdir}/seqkit.stats.tsv"
    """
}
