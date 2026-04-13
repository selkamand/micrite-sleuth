#!/usr/bin/env nextflow

// QC short alignments using picard and mosdepth
process QC_SHORT_ALIGNMENTS {
    tag "${meta.prefix}"

    input:
    // prefix should be "${sampleid}.${taxid}.${ref_id}"
    tuple val(meta), path(bam), path(bai)

    output:
    // tuple val(meta), path("*.mosdepth.global.dist.txt"), path("*.mosdepth.summary.txt"), emit: mosdepth
    // tuple val(meta), path("*.picard.alignment_summary_metrics.txt"), path("*.picard.insert_size_histogram.pdf"), path("*.picard.insert_size_metrics.txt"), emit: picard
    tuple val(meta), path("${meta.prefix}"), emit: multiqc
    tuple val(meta), path("${meta.prefix}/*"), emit: stats

    script:
    """
    set -euo pipefail

        outdir="${meta.prefix}"
        mkdir -p \${outdir}

# Picard: alignment summary
    picard CollectAlignmentSummaryMetrics \\
        I="${bam}" \\
        O="\${outdir}/${meta.prefix}.picard.alignment_summary_metrics.txt"

# Picard: insert size metrics (paired-end; will still run for SE but is less meaningful)
        picard CollectInsertSizeMetrics \\
        I="${bam}" \\
        O="\${outdir}/${meta.prefix}.picard.insert_size_metrics.txt" \\
        H="\${outdir}/${meta.prefix}.picard.insert_size_histogram.pdf"

# mosdepth: fast summaries only 
# (slowest part will be computing regions.bed.gz which computes coverage every 500bp window)
        mosdepth --by 500 -t ${task.cpus} -n "\${outdir}/${meta.prefix}" "${bam}"
# keeps: 
# \${outdir}/${meta.prefix}.mosdepth.summary.txt
# \${outdir}/${meta.prefix}.mosdepth.global.dist.txt
# \${outdir}/${meta.prefix}.mosdepth.regions.bed.gz
    """
}
