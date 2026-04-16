#!/usr/bin/env nextflow

// Use minimap2 to align de novo assembly to each reference genome
process ALIGN_WHOLE_GENOMES {

    // Preset for minimap alignment. See minimap help for -x argument for possible options
    ext preset: "asm5"

    input:
    tuple val(sampleid), val(taxid), path(assembly_fasta), val(ref_id), path(ref_fasta), path(ref_fai)

    output:
    tuple val(sampleid), val(taxid), val(ref_id), path(ref_fai), path("${sampleid}.${taxid}.${ref_id}.paf"), emit: full
    tuple val(sampleid), val(taxid), val(ref_id), path("${sampleid}.${taxid}.${ref_id}.paf"), emit: topublish
    tuple val(sampleid), val(taxid), path(ref_fasta), val(ref_id), path(assembly_fasta), path("${sampleid}.${taxid}.${ref_id}.paf"), emit: fordotplots

    script:
    """
    set -euo pipefail

    prefix="${sampleid}.${taxid}.${ref_id}"
    echo "Aligning assembly to ${ref_fasta} -> \${prefix}.paf" >&2

    minimap2 -t ${task.cpus} -cx ${task.ext.preset} --cs "${ref_fasta}" "${assembly_fasta}" > "\${prefix}.paf"
    """
}


process ALIGN_WHOLE_GENOMES_OLD {

    // Preset for minimap alignment. See minimap help for -x argument for possible options
    ext preset: "asm5"

    input:
    tuple val(sampleid), val(taxid), path(assembly_fasta), path(refgenomes)

    output:
    tuple val(sampleid), val(taxid), path("whole_genome_alignments")

    script:
    """
    set -euo pipefail

    outdir="whole_genome_alignments"
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

        echo "Aligning assembly to \$ref -> \${prefix}.paf" >&2

        minimap2 -t ${task.cpus} -cx ${task.ext.preset} --cs "\${ref}" "${assembly_fasta}" > "\${prefix}.paf"

    done

    # Optional: manifest file for easy downstream consumption
    ls -1 "\$outdir"/*.paf > "\$outdir/pafs.list"
    """
}
