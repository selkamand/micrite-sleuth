#!/usr/bin/env nexftlow

// Run FASTQC and seqkit stats on a pair of reads belonging to a particular taxid. 
// Assumes reads ahave been extracted by EXTRACT_READS_BY_TAXID

process QC_EXTRACTED_READS {
    tag "${sampleid}.${taxid}"

    input:
    tuple val(sampleid), val(taxid), path(fq1), path(fq2)

    output:
    tuple val(sampleid), val(taxid), path("seqkit.stats.tsv"), path("test.taxid_${taxid}.R1_fastqc.zip"), path("test.taxid_${taxid}.R2_fastqc.zip"), path("test.taxid_${taxid}.R1_fastqc.html"), path("test.taxid_${taxid}.R2_fastqc.html")

    script:
    """
    fastqc --nogroup ${fq1} ${fq2}
    seqkit stats ${fq1} ${fq2} > "seqkit.stats.tsv"
    """
}
