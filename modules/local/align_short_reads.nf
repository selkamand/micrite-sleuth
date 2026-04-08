#!/usr/bin/env nexftlow

// Use minimap2 to align short reads to a particular refgenome
process ALIGN_SHORT_READS_TO_GENOME {
    input:
    tuple val(sampleid), val(taxid), path(r1), path(r2), val(ref_id), path(ref_fasta), path(ref_fai)

    output:
    tuple val([sampleid: "${sampleid}", taxid: "${taxid}", ref_id: "${ref_id}", prefix: "${sampleid}.${taxid}.${ref_id}"]), path("${sampleid}.${taxid}.${ref_id}.sorted.bam"), path("${sampleid}.${taxid}.${ref_id}.sorted.bam.bai")

    script:
    """
    set -euo pipefail


    echo "Aligning reads to ${ref_fasta} -> ${sampleid}.${taxid}.${ref_id}.sorted.bam" >&2
    
    minimap2 -t ${task.cpus} -ax sr "${ref_fasta}" "${r1}" "${r2}" \\
    | samtools sort -@ ${task.cpus} -o "${sampleid}.${taxid}.${ref_id}.sorted.bam" -
        
    samtools index -@ ${task.cpus} "${sampleid}.${taxid}.${ref_id}.sorted.bam"

    """
}
