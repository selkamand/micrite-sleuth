#!/usr/bin/env nexftlow

// Assemble genome de novo using SPADES.py 
process ASSEMBLE {
    input:
    tuple val(sampleid), val(taxid), path(r1), path(r2)

    output:
    tuple val(sampleid), val(taxid), path("denovo_assembly_taxid_${taxid}/contigs.fasta"), path("denovo_assembly_taxid_${taxid}")

    script:
    """
    set -euo pipefail

    spades.py --threads ${task.cpus} --isolate -1 ${r1} -2 ${r2} -o "denovo_assembly_taxid_${taxid}"
    """
}
