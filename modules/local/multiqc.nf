process MULTIQC {
    tag "${meta.sampleid}.${meta.taxid}"

    container "community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b"

    input:
    tuple val(meta), path(alignment_stats)

    output:
    tuple val(meta), path("short_alignment_multiqc_report.html"), emit: report
    tuple val(meta), path("short_alignment_multiqc_report_data"), emit: data

    script:
    """
    set -euo pipefail

    multiqc \\
      --force \\
      --filename short_alignment_multiqc_report.html \\
      "${alignment_stats}"
    """
}

process MULTIQC_FILES {
    tag "${sampleid}.${taxid}"

    container "community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b"

    input:
    tuple val(sampleid), val(taxid), path(files)

    output:
    tuple val(sampleid), val(taxid), path("short_alignment_multiqc_report.html"), emit: report
    tuple val(sampleid), val(taxid), path("short_alignment_multiqc_report_data"), emit: data
    tuple val(sampleid), val(taxid), path("qcstats"), emit: stats

    script:
    """
    set -euo pipefail

    # Inputs arrive as multiple per-reference QC directories staged by Nextflow.
    # We copy all metrics into a single `qcstats/` dir so MultiQC runs on the combined set.
    mkdir -p qcstats

    shopt -s nullglob extglob
    # Copy everything from all top-level staged QC directories, excluding the output `qcstats/` itself.
    cp -p -- !(qcstats)/* qcstats/ || true

    multiqc \\
      --force \\
      --filename short_alignment_multiqc_report.html \\
      qcstats
    """
}
