process MULTIQC {
    tag "${sampleid}.${taxid}"

    container "community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b"

    input:
    tuple val(sampleid), val(taxid), path(alignment_stats)

    output:
    tuple val(sampleid), val(taxid), path("short_alignment_multiqc_report.html"), emit: report
    tuple val(sampleid), val(taxid), path("short_alignment_multiqc_report_data"), emit: data

    script:
    """
    set -euo pipefail

    multiqc \\
      --force \\
      --filename short_alignment_multiqc_report.html \\
      "${alignment_stats}"
    """
}
