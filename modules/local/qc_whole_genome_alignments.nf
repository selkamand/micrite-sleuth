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
process QC_WHOLE_GENOME_ALIGNMENTS_OLD {
    tag "${sampleid}.${taxid}"

    input:
    tuple val(sampleid), val(taxid), path(alndir)

    output:
    tuple val(sampleid), val(taxid), path("whole_genome_alignment_stats")

    script:
    """
    set -euo pipefail

    out="whole_genome_alignment_stats"
    mkdir -p "\$out"

    shopt -s nullglob

    pafs=( "${alndir}"/*.paf )

    if (( \${#pafs[@]} == 0 )); then
      echo "ERROR: No *paf files found in: ${alndir}" >&2
      exit 1
    fi

    for paf in "\${pafs[@]}"; do
      base=\$(basename "\$paf")
      prefix="\$out/\${base%.paf}"

      paftools.js stat "\${paf}" > \${prefix}.stats.tsv
      # paftools.js asmstat -q 0 -k 10000 -d 0.01 <ref.fai> "\${paf}" > \${prefix}.stats.tsv
    done

    """
}
