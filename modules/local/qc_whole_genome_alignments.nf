#!/usr/bin/env nexftlow

// use paftools to generate whole-genome alignment stats
process QC_WHOLE_GENOME_ALIGNMENTS {
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
