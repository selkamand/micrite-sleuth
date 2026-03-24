#!/usr/bin/env nexftlow

// Run FASTQC and seqkit stats on a pair of reads belonging to a particular taxid. 
// Assumes reads ahave been extracted by EXTRACT_READS_BY_TAXID

process QC_EXTRACTED_READS {
    tag "${sampleid}.${taxid}"
    cpus 2
    memory 512.MB

    input:
    tuple val(sampleid), val(taxid), path(fq1), path(fq2)

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.taxid_${taxid}.stats.tsv"), emit: seqkit
    tuple path("test.taxid_${taxid}.R1_fastqc.zip"), path("test.taxid_${taxid}.R2_fastqc.zip"), path("test.taxid_${taxid}.R1_fastqc.html"), path("test.taxid_${taxid}.R2_fastqc.html"), emit: fastqc

    script:
    """
    mkdir -p cache
    export XDG_CACHE_HOME=cache
    fastqc --memory 512MB -t ${task.cpus} --nogroup ${fq1} ${fq2}
    seqkit stats --threads ${task.cpus} --tabular ${fq1} ${fq2} > "${sampleid}.taxid_${taxid}.stats.tsv"
    """
}
