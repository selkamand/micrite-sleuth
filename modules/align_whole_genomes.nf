#!/usr/bin/env nexftlow

// Use minimap2 to align de novo assembly to each reference genome

process ALIGN_WHOLE_GENOMES {

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
