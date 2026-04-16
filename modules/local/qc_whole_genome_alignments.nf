#!/usr/bin/env nextflow



// use paftools to generate whole-genome alignment stats
process QC_WHOLE_GENOME_ALIGNMENTS {
    tag "${sampleid}.${taxid}"

    input:
    tuple val(sampleid), val(taxid), val(ref_id), path(fai), path(paf)

    output:
    tuple val(sampleid), val(taxid), path("*stats.tsv")

    script:
    """
    set -euo pipefail
    
      prefix="${sampleid}.${taxid}.${ref_id}"
      paftools.js stat "${paf}" > \${prefix}.stats.tsv
      paftools.js asmstat -q 0 -k 10000 -d 0.01 ${fai} ${paf} > \${prefix}.asmstats.tsv

    """
}


// Generate Dgenies index files 
process WGA_DOTPLOTS {
    tag "${sampleid}.${taxid}"
    container "selkamandcci/dgenies:0.0.1"

    input:
    tuple val(sampleid), val(taxid), path(ref_fasta), val(ref_id), path(assembly_fasta), path(paf)

    output:
    tuple val(sampleid), val(taxid), path("${ref_id}.dgenies.idx"), path("${sampleid}.${taxid}.denovo_assembly.dgenies.idx")

    script:
    """
    set -euo pipefail
   
    # Index Reference Fasta 
    index_fasta.py -i "${ref_fasta}" -n ${ref_id} -o ${ref_id}.dgenies.idx
    
    # Index Query Fasta
    index_fasta.py -i "${assembly_fasta}" -n "${sampleid}.${taxid}" -o "${sampleid}.${taxid}.denovo_assembly.dgenies.idx"
    
    # Run Dgenies on PAF 
    # TODO: run standalone dgenies locally
    """
}
