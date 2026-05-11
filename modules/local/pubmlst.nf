process PUBMLST {
    tag "${sampleid}|${taxid}"

    input:
    tuple val(sampleid), val(taxid), path(assembly_fasta)
    path database

    output:
    tuple val(sampleid), val(taxid), path("pubmlst.csv")

    script:
    """
    set -euo pipefail

    mlst --db ${database} "${assembly_fasta}" > pubmlst.csv
    """
}
