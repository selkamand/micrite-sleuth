#!/usr/bin/env nexftlow

// Use bowtie2 to align short reads to a particular refgenome
process ALIGN_SHORT_READS_TO_GENOME {
    input:
    tuple val(sampleid), val(taxid), path(r1), path(r2), path(refgenomes)

    output:
    tuple val(sampleid), val(taxid), path("short_read_alignments")

    script:
    """
    set -euo pipefail

    outdir="short_read_alignments"
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
