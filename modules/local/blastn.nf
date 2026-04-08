#!/usr/bin/env nextflow

// Use fq to substample nreads then convert to fasta to run blast
process BLASTN {

    container "community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759"
    cpus 1

    input:
    tuple val(sampleid), val(taxid), path(fasta)

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.taxid_${taxid}.blastn.tsv")

    script:
    """
    set -euo pipefail
    blastn -db nt \
    -evalue 1e-10 \
    -max_target_seqs 20 \
    -outfmt '7 std staxid qcovs qcovhsp stitle' \
    -perc_identity 90 -max_hsps 10 -subject_besthit \
    -remote \
    -out ${sampleid}.taxid_${taxid}.blastn.tsv \
    -query  ${fasta}
    """
}

process BLASTPARSE_RUN {

    container "selkamandcci/blastparse:rocker_0.0.1"

    cpus 1

    ext evalue: "1e-10", pident: 90

    input:
    tuple val(sampleid), val(taxid), path(fasta)

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.taxid_${taxid}.blastn.tsv"), path("${sampleid}.taxid_${taxid}.blastn.config.tsv")

    script:
    """
    set -euo pipefail

    Rscript -e 'blastparse::blast_run(
        query = "${fasta}", 
        db = "nt", 
        remote = TRUE,
        overwrite = TRUE,
        evalue = "${task.ext.evalue}",
        perc_identity = ${task.ext.pident},
        outfile_prefix = "${sampleid}.taxid_${taxid}"
    )'
    """
}
