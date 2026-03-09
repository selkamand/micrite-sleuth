#!/usr/bin/env nexftlow

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

    // output directory
    outdir: Path = "micritesleuth"
}

include { EXTRACT_READS_BY_TAXID } from "./modules/krakentools.nf"
include { QC_EXTRACTED_READS } from "./modules/qc_extracted_reads.nf"
include { ALIGN_SHORT_READS_TO_GENOME } from "./modules/align_short_reads.nf"
include { QC_SHORT_ALIGNMENTS } from "./modules/qc_short_alignments.nf"

workflow {

    main:
    def r1 = file(params.r1)
    def r2 = file(params.r2)
    def kraken = file(params.kraken)
    def kreport = file(params.kreport)
    def sampleid = params.sampleid
    def refgenomes = file(params.refgenomes)

    // Extract reads that hit taxid (or descendant) 
    reads_from_taxid_ch = EXTRACT_READS_BY_TAXID(channel.of(tuple(sampleid, params.taxid, kraken, kreport, r1, r2)))

    // QC the extracted reads
    qc_from_taxid_ch = QC_EXTRACTED_READS(reads_from_taxid_ch)

    // Align extracted reads to each reference genome in the provided directory
    align_in_ch = reads_from_taxid_ch.map { sid, tx, fq1, fq2 ->
        tuple(sid, tx, fq1, fq2, refgenomes)
    }
    aligned_to_refs_ch = ALIGN_SHORT_READS_TO_GENOME(align_in_ch)

    // QC the short read alignments with picard and mosdepth
    qc_short_alignments_ch = QC_SHORT_ALIGNMENTS(aligned_to_refs_ch)

    publish:
    extracted_reads = reads_from_taxid_ch
    extracted_read_fastqc = qc_from_taxid_ch
    alignments = aligned_to_refs_ch
    alignment_stats = qc_short_alignments_ch
}

output {
    extracted_reads {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/reads/"
        mode 'copy'
    }
    extracted_read_fastqc {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/read_stats"
        mode 'copy'
    }
    alignments {
        path "${params.outdir}/${params.sampleid}/${params.taxid}"
        mode 'copy'
    }
    alignment_stats {
        path "${params.outdir}/${params.sampleid}/${params.taxid}"
        mode 'copy'
    }
}
