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

    // Size of sleuthing
    outdir: Path = "micritesleuth"
}

process EXTRACT_READS_BY_TAXIDS {
    input:
    tuple val(sampleid), val(taxid), path(kraken), path(kreport), path(r1), path(r2)

    output:
    tuple val(sampleid), path("${sampleid}.${taxid}.R1.fq"), path("${sampleid}.${taxid}.R2.fq")

    script:
    """
    extract_kraken_reads.py \
    --include-children -k ${kraken} \
    --fastq-output \
    -r ${kreport} \
    -t ${taxid} \
    -s ${r1} -s2 ${r2} \
    -o ${sampleid}.${taxid}.R1.fq -o2 ${sampleid}.${taxid}.R2.fq 
  """
}

workflow {

    main:
    def r1 = file(params.r1)
    def r2 = file(params.r2)
    def kraken = file(params.kraken)
    def kreport = file(params.kreport)
    def sampleid = params.sampleid

    reads_from_taxid_ch = EXTRACT_READS_BY_TAXIDS(channel.of(tuple(sampleid, params.taxid, kraken, kreport, r1, r2)))

    publish:
    extracted_reads = EXTRACT_READS_BY_TAXIDS.out
}

output {
    extracted_reads {
        path "${params.outdir}/${params.sampleid}/"
        mode 'copy'
    }
}
