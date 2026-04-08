#!/usr/bin/env nexftlow

// QC short alignments using picard and mosdepth
process QC_SHORT_ALIGNMENTS_OLD {
    tag "${sampleid}.${taxid}"

    input:
    tuple val(sampleid), val(taxid), path(alndir)

    output:
    tuple val(sampleid), val(taxid), path("alignment_stats")

    script:
    """
    set -euo pipefail

    out="alignment_stats"
    mkdir -p "\$out"

    shopt -s nullglob

    bams=( "${alndir}"/*.sorted.bam )

    if (( \${#bams[@]} == 0 )); then
      echo "ERROR: No *.sorted.bam files found in: ${alndir}" >&2
      exit 1
    fi

    for bam in "\${bams[@]}"; do
      base=\$(basename "\$bam")
      prefix="\$out/\${base%.sorted.bam}"

      # Basic sanity: ensure index exists (samtools index creates .bam.bai by default)
      if [[ ! -e "\$bam.bai" ]]; then
        echo "ERROR: Missing BAM index: \$bam.bai" >&2
        exit 1
      fi

      # Picard: alignment summary
      picard CollectAlignmentSummaryMetrics \\
        I="\$bam" \\
        O="\${prefix}.picard.alignment_summary_metrics.txt"

      # Picard: insert size metrics (paired-end; will still run for SE but is less meaningful)
      picard CollectInsertSizeMetrics \\
        I="\$bam" \\
        O="\${prefix}.picard.insert_size_metrics.txt" \\
        H="\${prefix}.picard.insert_size_histogram.pdf"

      # mosdepth: fast summaries only
      mosdepth -t ${task.cpus} -n "\${prefix}.mosdepth" "\$bam"
      # keeps: \${prefix}.mosdepth.summary.txt and \${prefix}.mosdepth.global.dist.txt (and a few small extras)
    done

    """
}


process QC_SHORT_ALIGNMENTS {
    tag "${meta.prefix}"

    input:
    // prefix should be "${sampleid}.${taxid}.${ref_id}"
    tuple val(meta), path(bam), path(bai)

    output:
    // tuple val(meta), path("*.mosdepth.global.dist.txt"), path("*.mosdepth.summary.txt"), emit: mosdepth
    // tuple val(meta), path("*.picard.alignment_summary_metrics.txt"), path("*.picard.insert_size_histogram.pdf"), path("*.picard.insert_size_metrics.txt"), emit: picard
    tuple val(meta), path("qc"), emit: multiqc

    script:
    """
    set -euo pipefail

        outdir="qc"
        mkdir -p \${outdir}

# Picard: alignment summary
    picard CollectAlignmentSummaryMetrics \\
        I="${bam}" \\
        O="\${outdir}.${meta.prefix}.picard.alignment_summary_metrics.txt"

# Picard: insert size metrics (paired-end; will still run for SE but is less meaningful)
        picard CollectInsertSizeMetrics \\
        I="${bam}" \\
        O="\${outdir}.${meta.prefix}.picard.insert_size_metrics.txt" \\
        H="\${outdir}.${meta.prefix}.picard.insert_size_histogram.pdf"

# mosdepth: fast summaries only
        mosdepth -t ${task.cpus} -n "\${outdir}.${meta.prefix}" "${bam}"
# keeps: \${outdir}.${meta.prefix}.mosdepth.summary.txt and \${outdir}.${meta.prefix}.mosdepth.global.dist.txt (and a few small extras)

    """
}
