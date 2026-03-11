#!/usr/bin/env nextflow

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

    // [[ Conditional parts of pipelines ]]

    // Run remote BLAST
    run_remote_blastn: Boolean = false

    // Number of reads to remote blast 
    blastn_reads: Integer = 200

    // Run de novo assembly
    run_assembly: Boolean = true

    // Run bacterial genome annotation
    run_assembly_annotation: Boolean = false

    // Path to bakta database directory
    bakta_database: Path? = null

    // QUAST de novo genome QC config 
    quast_min_contig: Integer
    quast_est_ref_size: Integer


    // output directory
    outdir: Path = "micritesleuth"
}

include { EXTRACT_READS_BY_TAXID } from "./modules/krakentools.nf"
include { QC_EXTRACTED_READS } from "./modules/qc_extracted_reads.nf"
include { ALIGN_SHORT_READS_TO_GENOME } from "./modules/align_short_reads.nf"
include { QC_SHORT_ALIGNMENTS } from "./modules/qc_short_alignments.nf"
include { SUBSAMPLE_FQ } from "./modules/subsample.nf"
include { BLASTN } from "./modules/blastn.nf"
include { ASSEMBLE } from "./modules/assemble.nf"
include { ALIGN_WHOLE_GENOMES } from "./modules/align_whole_genomes.nf"
include { QC_WHOLE_GENOME_ALIGNMENTS } from "./modules/qc_whole_genome_alignments.nf"
include { ANNOTATE_BACTERIAL_GENOME } from "./modules/annotate_bacterial_genome.nf"
include { BARRNAP } from "./modules/extract_rrna_seqs.nf"
include { QUAST_WHOLE_GENOME_ASSEMBLY } from "./modules/quast.nf"

workflow {

    main:
    def r1 = file(params.r1)
    def r2 = file(params.r2)
    def kraken = file(params.kraken)
    def kreport = file(params.kreport)
    def sampleid = params.sampleid
    def refgenomes = file(params.refgenomes)
    def bakta_database = params.bakta_database != null ? file(params.bakta_database) : null

    //Parameter checks
    if (params.run_assembly_annotation) {
        if (bakta_database == null) {
            error("--bakta_database paramater must be provided when --run_assembly_annotation is true")
        }
        if (!bakta_database.exists()) {
            error("Failed to find bakta genome annotation database: [${params.bakta_database}]. Please ensure --bakta_database parameter points to a valid bakta database.")
        }
    }


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

    // TODO: update this so we dynamically check that taxid is indeed bacterial 
    taxid_is_bacterial = true

    // Initialize optional outputs as empty channels
    subsample_for_blastn_ch = channel.empty()
    blastn_ch = channel.empty()
    assembly_ch = channel.empty()
    annotation_ch = channel.empty()
    quast_ch = channel.empty()

    // Subsample and remote blast 
    if (params.run_remote_blastn) {

        // Define how many reads to blastn (will be taken from R1 fq)
        def nreads = params.blastn_reads
        // TODO: if user-specified blastn_nreads is > than the number of reads in R1 fastq, we should just blast all possible reads in FQ.

        // Perform subsampling
        subsample_for_blastn_ch = reads_from_taxid_ch.map { sid, tx, fq1, _fq2 -> tuple(sid, tx, fq1, nreads) }
            | SUBSAMPLE_FQ

        // Perform blastn 
        blastn_ch = BLASTN(subsample_for_blastn_ch)
    }

    // Run de novo assembly 
    if (params.run_assembly) {
        // TODO: add subsample to 30x depth if we have enough reads 

        // Create de novo assembly
        // Note assembly_raw channel includes both contigs.fasta AND the whole genome Dir
        assembly_raw_ch = ASSEMBLE(reads_from_taxid_ch)

        // Create assembly output channel (we don't want to output assembly fasta twice, once from direct path & once from directory)
        assembly_output_ch = assembly_ch.map { sid, tx, _assembly_fasta, assembly_dir -> tuple(sid, tx, assembly_dir) }

        // Create a simpler channel that just includes sample, taxid, and path to fasta
        // This is the channel we'll feed into most downstream operations
        assembly_ch = assembly_raw_ch.map { sid, tx, assembly_fasta, _assembly_dir -> tuple(sid, tx, assembly_fasta) }

        // Perform whole-genome alignments
        whole_genome_alignments_ch = assembly_ch.map { sid, tx, assembly_fasta -> tuple(sid, tx, assembly_fasta, refgenomes) }
            | ALIGN_WHOLE_GENOMES

        // Compute Stats on whole genome alignments
        whole_genome_alignment_stats_ch = QC_WHOLE_GENOME_ALIGNMENTS(whole_genome_alignments_ch)

        // Barrnap 16/23S rRNA extraction from de ovo assembly
        barrnap_ch = BARRNAP(assembly_ch)

        // Run QUAST QC on de novo assembly.
        quast_in = assembly_ch.map { sid, tx, assembly_fasta ->
            tuple(
                sid,
                tx,
                assembly_fasta,
                params.quast_min_contig,
                params.quast_est_ref_size,
            )
        }

        quast_ch = QUAST_WHOLE_GENOME_ASSEMBLY(quast_in)

        // Annotate genome with BAKTA
        if (params.run_assembly_annotation & taxid_is_bacterial) {
            annotation_ch = assembly_ch.map { sid, tx, assembly_fasta -> tuple(sid, tx, assembly_fasta, bakta_database) }
                | ANNOTATE_BACTERIAL_GENOME
        }
    }

    publish:
    extracted_reads = reads_from_taxid_ch
    extracted_read_fastqc = qc_from_taxid_ch
    alignments = aligned_to_refs_ch
    alignment_stats = qc_short_alignments_ch
    subsampled_reads_for_blastn = subsample_for_blastn_ch
    blastn = blastn_ch
    assembly = assembly_output_ch
    whole_genome_alignments = whole_genome_alignments_ch
    whole_genome_alignment_stats = whole_genome_alignment_stats_ch
    annotation = annotation_ch
    barrnap = barrnap_ch
    quast = quast_ch
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
    subsampled_reads_for_blastn {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/blastn/"
        mode 'copy'
    }
    blastn {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/blastn/"
        mode 'copy'
    }
    assembly {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/"
        mode 'copy'
    }
    whole_genome_alignments {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/"
        mode 'copy'
    }
    whole_genome_alignment_stats {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/"
        mode 'copy'
    }
    barrnap {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/barrnap/"
        mode 'copy'
    }
    annotation {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/"
        mode 'copy'
    }
    quast {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/"
        mode 'copy'
    }
}
