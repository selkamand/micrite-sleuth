#!/usr/bin/env nexftlow

// Use fq to substample nreads then convert to fasta to run blast
process BLASTN {
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
    -outfmt '6 std staxid qcovs qcovhsp stitle' \
    -perc_identity 90 -max_hsps 10 -subject_besthit \
    -remote \
    -out ${sampleid}.taxid_${taxid}.blastn.tsv \
    -query  ${fasta}
    """
}
