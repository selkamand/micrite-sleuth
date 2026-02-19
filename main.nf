params {

    // sample identifier
    sampleid: String

    //Kraken result file (stdout)
    kraken: Path

    // Kraken Report (used to identify child taxids)
    kreport: Path

    // Path to fastq reads classified
    r1: Path

    // Path to fastq reads classified
    r2: Path

    // Taxid to extract (will include children)
    taxid: Integer

    // A folder containing bwa-mem2 + minimap2 indexed microbial genomes 
    // Reads assigned to taxid X will be aligned to those reference genomes
    refgenomes: Path


    // Size of sleuthing
    outdir: Path = "micritesleuth"
}

process EXTRACT_READS_BY_TAXIDS {
    input:
    tuple val(sampleid), val(taxid), path(kraken), path(kreport), path(r1), path(r2)

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.taxid_${taxid}.R1.fq"), path("${sampleid}.taxid_${taxid}.R2.fq")

    script:
    """
    extract_kraken_reads.py \
    --include-children -k ${kraken} \
    --fastq-output \
    -r ${kreport} \
    -t ${taxid} \
    -s ${r1} -s2 ${r2} \
    -o ${sampleid}.taxid_${taxid}.R1.fq -o2 ${sampleid}.taxid_${taxid}.R2.fq 
  """
}


process QC_EXTRACTED_READS {
    tag "${sampleid}.${taxid}"

    input:
    tuple val(sampleid), val(taxid), path(fq1), path(fq2)

    output:
    tuple val(sampleid), val(taxid), path("stats")

    script:
    """
    mkdir -p stats
    fastqc --nogroup -o stats ${fq1} ${fq2}
    seqkit stats ${fq1} ${fq2} > stats/seqkit.stats.tsv
    """
}

process ALIGN_SHORT_READS_TO_GENOME {
    input:
    tuple val(sampleid), val(taxid), path(r1), path(r2), path(refgenomes)

    output:
    tuple val(sampleid), val(taxid), path("aligning_${taxid}_reads_to_refgenomes")

    script:
    """
    set -euo pipefail

    outdir="aligning_${taxid}_reads_to_refgenomes"
    mkdir -p "\$outdir"

    shopt -s nullglob

    # Collect references (support common FASTA extensions; include gz)
    refs=( "${refgenomes}"/*.fa "${refgenomes}"/*.fna "${refgenomes}"/*.fasta "${refgenomes}"/*.fa.gz "${refgenomes}"/*.fna.gz "${refgenomes}"/*.fasta.gz )

    if (( \${#refs[@]} == 0 )); then
        echo "ERROR: No reference genomes found in: ${refgenomes}" >&2
        echo "Looked for: *.fa *.fna *.fasta (optionally .gz)" >&2
        exit 1
    fi

    echo "Found \${#refs[@]} reference genome(s) in ${refgenomes}" >&2

    for ref in "\${refs[@]}"; do
        base=\$(basename "\$ref")
        # strip extensions safely
        name="\${base%.gz}"
        name="\${name%.fasta}"
        name="\${name%.fna}"
        name="\${name%.fa}"

        prefix="\$outdir/${sampleid}.${taxid}.\$name"

        echo "Aligning reads to \$ref -> \${prefix}.sorted.bam" >&2

        # minimap2 outputs SAM to stdout; samtools sort makes coordinate-sorted BAM
        minimap2 -t ${task.cpus} -ax sr "\$ref" "${r1}" "${r2}" \\
          | samtools sort -@ ${task.cpus} -o "\${prefix}.sorted.bam" -

        samtools index -@ ${task.cpus} "\${prefix}.sorted.bam"

        
    done

    # Optional: manifest file for easy downstream consumption
    ls -1 "\$outdir"/*.sorted.bam > "\$outdir/bams.list"
    """
}
process SHORT_ALIGNMENT_STATS {
    tag "${sampleid}.${taxid}"

    input:
    tuple val(sampleid), val(taxid), path(alndir)

    output:
    tuple val(sampleid), val(taxid), path("alignment_stats_${taxid}")

    script:
    """
    set -euo pipefail

    out="alignment_stats_${taxid}"
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

workflow {

    main:
    def r1 = file(params.r1)
    def r2 = file(params.r2)
    def kraken = file(params.kraken)
    def kreport = file(params.kreport)
    def sampleid = params.sampleid
    def refgenomes = file(params.refgenomes)

    // Extract reads that hit taxid (or descendant) 
    reads_from_taxid_ch = EXTRACT_READS_BY_TAXIDS(channel.of(tuple(sampleid, params.taxid, kraken, kreport, r1, r2)))

    // QC the extracted reads
    qc_from_taxid_ch = QC_EXTRACTED_READS(reads_from_taxid_ch)

    // Align extracted reads to each reference genome in the provided directory
    align_in_ch = reads_from_taxid_ch.map { sid, tx, fq1, fq2 ->
        tuple(sid, tx, fq1, fq2, refgenomes)
    }
    aligned_to_refs_ch = ALIGN_SHORT_READS_TO_GENOME(align_in_ch)
    alignment_stats_ch = SHORT_ALIGNMENT_STATS(aligned_to_refs_ch)

    publish:
    extracted_reads = EXTRACT_READS_BY_TAXIDS.out
    extracted_fastqc = qc_from_taxid_ch
    alignments = aligned_to_refs_ch
    alignment_stats = alignment_stats_ch
}

output {
    extracted_reads {
        path "${params.outdir}/${params.sampleid}/"
        mode 'copy'
    }
    extracted_fastqc {
        path "${params.outdir}/${params.sampleid}/"
        mode 'copy'
    }
    alignments {
        path "${params.outdir}/${params.sampleid}/"
        mode 'copy'
    }
    alignment_stats {
        path "${params.outdir}/${params.sampleid}/"
        mode 'copy'
    }
}
