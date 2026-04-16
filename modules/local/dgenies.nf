// Create denies index for fasta file
process PREPARE_DGENIES_INDEX {

    tag "${ref_id}"
    container "selkamandcci/dgenies:0.0.1"

    input:
    tuple val(ref_id), path(ref_fasta), path(ref_fai)

    output:
    tuple val(ref_id), path("${ref_id}.dgenies.idx")

    script:
    """
    set -euo pipefail

    index_fasta.py -i "${ref_fasta}" -n ${ref_id} -o ${ref_id}.dgenies.idx
    """
}
