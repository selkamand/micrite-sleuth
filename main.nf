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
    skip_pre_assembly_downsampling: Boolean = false

    // Run bacterial genome annotation
    run_assembly_annotation: Boolean = false

    // Path to bakta database directory
    bakta_database: Path? = null

    // What size would we expect the genome size to be 
    // (used in quast QC NG50 calculations and subsampling reads to ~30X coverage ahead of assembly)
    genome_size_guess: Integer

    // QUAST de novo genome QC config 
    quast_min_contig: Integer

    // Busco settings for genome completeness
    run_busco: Boolean = false
    busco_lineage: String?
    busco_dataset: Path?

    // Silva ribosomal RNA databases. Used to classify any 16S/23S rRNA sequences extracted by barrnap
    arb_16s: Path
    arb_23s: Path

    // output directory
    outdir: Path = "micritesleuth"
}

include { EXTRACT_READS_BY_TAXID } from "./modules/local/krakentools.nf"
include { QC_EXTRACTED_READS } from "./modules/local/qc_extracted_reads.nf"
include { ALIGN_SHORT_READS_TO_GENOME } from "./modules/local/align_short_reads.nf"
include { QC_SHORT_ALIGNMENTS } from "./modules/local/qc_short_alignments.nf"
include { SUBSAMPLE_FQ } from "./modules/local/subsample.nf"
include { BLASTN } from "./modules/local/blastn.nf"
include { ASSEMBLE } from "./modules/local/assemble.nf"
include { ALIGN_WHOLE_GENOMES } from "./modules/local/align_whole_genomes.nf"
include { WGA_DOTPLOTS ; QC_WHOLE_GENOME_ALIGNMENTS } from "./modules/local/qc_whole_genome_alignments.nf"
include { ANNOTATE_BACTERIAL_GENOME } from "./modules/local/annotate_bacterial_genome.nf"
include { BARRNAP ; SPLIT_RRNA_FASTA } from "./modules/local/extract_rrna_seqs.nf"
include { SINA_SEARCH_AND_CLASSIFY } from "./modules/local/classify_rrna_seqs.nf"
include { QUAST_WHOLE_GENOME_ASSEMBLY } from "./modules/local/quast.nf"
include { COUNT_TOTAL_BASES } from "./modules/local/parse_seqkit_stats.nf"
include { SUBSAMPLE_BY_PROPORTION } from "./modules/local/subsample_by_proportion.nf"
include { BUSCO_COMPLETENESS } from "./modules/local/busco.nf"
include { MULTIQC_FILES } from './modules/local/multiqc.nf'
include { BLASTPARSE_RUN } from './modules/local/blastn.nf'


workflow {

    main:
    def r1 = file(params.r1)
    def r2 = file(params.r2)
    def kraken = file(params.kraken)
    def kreport = file(params.kreport)
    def sampleid = params.sampleid
    def refgenomes = file(params.refgenomes)
    def bakta_database = params.bakta_database != null ? file(params.bakta_database) : null
    def arb_16s = file(params.arb_16s)
    def arb_23s = file(params.arb_23s)

    // When downsampling reads for de novo assembly we should aim for ~30x coverage
    def assembly_target_cov = 30

    // Identify all reference genomes we want to compare with based on 
    def ref_glob = "${refgenomes}/*.{fa,fna,fasta,fa.gz,fna.gz,fasta.gz}"
    def ref_paths_ch = channel.fromPath(ref_glob, checkIfExists: true)

    // Define a channel of 3 element tuples. 
    // 0=reference genome id (from filename), 
    // 1=reference genome fasta file, 
    // 2=fai index
    def refgenomes_ch = ref_paths_ch.map { ref_fasta_path ->
        def fai = file("${ref_fasta_path}.fai")
        def ref_name = ref_fasta_path.getFileName().toString()
        def ref_id = ref_name
            .replaceFirst(/\.gz$/, '')
            .replaceFirst(/\.fasta$/, '')
            .replaceFirst(/\.fna$/, '')
            .replaceFirst(/\.fa$/, '')

        if (!fai.exists()) {
            error("Reference genome missing fai index. Please run: `samtools faidx ${ref_fasta_path}` then try again")
        }
        tuple(ref_id, ref_fasta_path, fai)
    }
    //.view()

    // if (ref_entries.isEmpty()) {
    // error("No reference FASTA files found at --refgenomes path using glob: ${ref_glob}")
    //}


    //Parameter checks
    if (params.run_assembly == false) {
        if (params.run_assembly_annotation) {
            error("Can NOT run assembly annotation if assembly is going to be skipped. Either set --run_assembly_annotation to false or --run_assembly to true")
        }
        if (params.run_busco) {
            error("Can NOT run assembly completeness (BUSCO) evaluation if assembly is going to be skipped. Either set --run_busco to false or --run_assembly to true")
        }
    }

    if (params.run_busco) {
        if (params.busco_lineage == null) {
            error("A --busco_lineage must be specified if --run_busco flag is supplied")
        }
        if (params.busco_dataset == null) {
            error("A busco dataset must be specified if --run_busco flag is supplied")
        }

        busco_dataset = file(params.busco_dataset)
        if (!busco_dataset.exists()) {
            error("BUSCO dataset ${busco_dataset} does not exist")
        }
    }

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

    // Create a pairwise combination of every set of extracted reads and every reference genome 
    // so we can align to each refgenome
    aligned_input_ch = reads_from_taxid_ch.combine(refgenomes_ch)

    // Align short reads to reference genomes
    aligned_to_refs_ch = ALIGN_SHORT_READS_TO_GENOME(aligned_input_ch)
    // def aligned_to_refs_ch = channel.empty()

    // QC the short read alignments with picard and mosdepth (runs once per ref genome)
    qc_short_alignments_ch = QC_SHORT_ALIGNMENTS(aligned_to_refs_ch)

    // MultiQC should run once per taxid/sample, combining QC outputs from all reference genomes.
    qc_short_alignments_grouped_for_multiqc_ch = qc_short_alignments_ch.multiqc
        .map { meta, qc_dir -> tuple(meta.sampleid, meta.taxid, qc_dir) }
        .groupTuple(by: [0, 1])

    // Compile short read alignment qc with multiqc (runs once per taxid/sample)
    multiqc_short_alignments_ch = MULTIQC_FILES(
        qc_short_alignments_grouped_for_multiqc_ch.map { sid, tx, qc_dirs ->
            tuple(sid, tx, qc_dirs)
        }
    )

    // TODO: update this so we dynamically check that taxid is indeed bacterial 
    taxid_is_bacterial = true

    // Initialize optional outputs as empty channels
    subsample_for_blastn_ch = channel.empty()
    blastn_ch = channel.empty()
    assembly_ch = channel.empty()
    annotation_ch = channel.empty()
    quast_ch = channel.empty()

    // Subsample a small number of reads classified at/under taxid  
    // (these will later be used for blastn)
    // Define how many reads to blastn (will be taken from R1 fq)
    // TODO: if user-specified blastn_nreads is > than the number of reads in R1 fastq, we should just blast all possible reads in FQ.
    def nreads = params.blastn_reads

    // Perform subsampling
    subsample_for_blastn_ch = reads_from_taxid_ch.map { sid, tx, fq1, _fq2 -> tuple(sid, tx, fq1, nreads) }
        | SUBSAMPLE_FQ

    // Perform blastn
    if (params.run_remote_blastn) {
        // blastn_ch = BLASTN(subsample_for_blastn_ch)
        blastn_ch = BLASTPARSE_RUN(subsample_for_blastn_ch)
    }
    else {
        blastn_ch = channel.empty()
    }

    // Run de novo assembly 
    if (params.run_assembly) {

        // Subsample reads to ~30x coverage (or whatever assembly_target_cov is set to) based on genome_size_guess. 
        // Parse total reads into nextflow channel (used later for downsampling to 30x estimated genome size)
        stats_for_downsampling_raw_ch = COUNT_TOTAL_BASES(qc_from_taxid_ch.seqkit).map { id, taxid, statfile ->
            def total_length = statfile.text.trim() as Integer
            def current_depth = total_length / params.genome_size_guess
            def subsample_prop = assembly_target_cov / current_depth
            // def stats = [total_length: total_length, current_depth: current_depth, subsample_prop: subsample_prop, target_cov: assembly_target_cov]
            tuple(id, taxid, subsample_prop)
        }
        //.v}iew { v -> "Extracted Read Stats for Downsampling: ${v}" }

        // Enrich with actual read paths 
        // Also branch based on wether we have too many / too few reads
        stats_for_downsampling_ch = reads_from_taxid_ch
            .join(stats_for_downsampling_raw_ch, by: [0, 1])
            .branch { _sid, _tx, _fq1, _fq2, prop ->
                subsample: prop < 1 && !params.skip_pre_assembly_downsampling
                full: true
            }

        // stats_for_downsampling_ch.subsample.view { v -> "[Subsample] ${v}" }
        // stats_for_downsampling_ch.full.view { v -> "[Branch] ${v}" }

        // Subsample reads (will be empty if branch is full instead of subsample)
        subsampled_reads_for_assembly_ch = SUBSAMPLE_BY_PROPORTION(stats_for_downsampling_ch.subsample)

        // If branch is full just pass along reads
        full_taxid_reads_for_assembly_ch = stats_for_downsampling_ch.full.map { sid, tx, fq1, fq2, _prop -> tuple(sid, tx, fq1, fq2) }

        // Create assembly imput by mixing the two possible channel branches
        assembly_input_ch = full_taxid_reads_for_assembly_ch.mix(subsampled_reads_for_assembly_ch)

        // Create de novo assembly
        // Note assembly_raw channel includes both contigs.fasta AND the whole genome Dir
        assembly_ch = ASSEMBLE(assembly_input_ch)

        // Perform whole-genome alignments against every refgenome
        whole_genome_alignments_ch = assembly_ch.contigs.combine(refgenomes_ch)
            | ALIGN_WHOLE_GENOMES

        // Compute Stats on whole genome alignments
        whole_genome_alignment_stats_ch = QC_WHOLE_GENOME_ALIGNMENTS(whole_genome_alignments_ch.full)


        // Run Dgenies on de novo assembly whole-genome alignments
        ch_dgenies = WGA_DOTPLOTS(whole_genome_alignments_ch.fordotplots)

        // Barrnap 16/23S rRNA extraction from de novo assembly
        barrnap_ch = BARRNAP(assembly_ch.contigs)

        // Split barrnap rRNAs to individual fastqs for small and large subunut components 
        rrna_sequences_ch = SPLIT_RRNA_FASTA(barrnap_ch)

        // Classify 16S sequence against SINA database
        sina_classified_16s_ch = SINA_SEARCH_AND_CLASSIFY(rrna_sequences_ch.SSU_16S, arb_16s, "16S")
        sina_classified_23s_ch = SINA_SEARCH_AND_CLASSIFY(rrna_sequences_ch.LSU_23S, arb_23s, "23S")

        // Run QUAST QC on de novo assembly.
        quast_in = assembly_ch.contigs.map { sid, tx, assembly_fasta ->
            tuple(
                sid,
                tx,
                assembly_fasta,
                params.quast_min_contig,
                params.genome_size_guess,
            )
        }

        quast_ch = QUAST_WHOLE_GENOME_ASSEMBLY(quast_in)

        // Annotate genome with BAKTA
        if (params.run_assembly_annotation & taxid_is_bacterial) {
            annotation_ch = assembly_ch.contigs.map { sid, tx, assembly_fasta -> tuple(sid, tx, assembly_fasta, bakta_database) }
                | ANNOTATE_BACTERIAL_GENOME
        }
        else {
            annotation_ch = channel.empty()
        }

        // BUSCO completeness
        if (params.run_busco) {
            busco_ch = BUSCO_COMPLETENESS(assembly_ch.contigs, params.busco_lineage, busco_dataset)
        }
        else {
            busco_ch = channel.empty()
        }
    }

    publish:
    extracted_reads = reads_from_taxid_ch
    extracted_read_fastqc = qc_from_taxid_ch.fastqc
    extracted_read_seqkit = qc_from_taxid_ch.seqkit
    short_read_alignments = aligned_to_refs_ch
    short_read_alignment_stats = multiqc_short_alignments_ch.stats
    short_read_alignment_multiqc_report = multiqc_short_alignments_ch.report
    short_read_alignment_multiqc_data = multiqc_short_alignments_ch.data
    subsampled_reads_for_blastn = subsample_for_blastn_ch
    blastn = blastn_ch
    subsampled_reads_for_assembly = subsampled_reads_for_assembly_ch
    assembly = assembly_ch.all_results
    whole_genome_alignments = whole_genome_alignments_ch.topublish
    whole_genome_alignment_stats = whole_genome_alignment_stats_ch
    annotation = annotation_ch
    barrnap = barrnap_ch
    ssu_16s = rrna_sequences_ch.SSU_16S
    lsu_23s = rrna_sequences_ch.LSU_23S
    lsu_5s = rrna_sequences_ch.LSU_5S
    sina_classified_16s = sina_classified_16s_ch.all
    sina_classified_23s = sina_classified_23s_ch.all
    quast = quast_ch
    busco = busco_ch
    dgenies = ch_dgenies
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
    extracted_read_seqkit {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/read_stats"
        mode 'copy'
    }
    short_read_alignments {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/short_read_alignments"
        mode 'copy'
    }
    short_read_alignment_stats {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/short_read_alignments/"
        mode 'copy'
    }
    short_read_alignment_multiqc_report {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/short_read_alignments/multiqc"
        mode 'copy'
    }
    short_read_alignment_multiqc_data {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/short_read_alignments/multiqc"
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
    subsampled_reads_for_assembly {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/reads/subsampled_for_assembly/"
        mode 'copy'
    }
    assembly {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/"
        mode 'copy'
    }
    whole_genome_alignments {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/whole_genome_alignments/"
        mode 'copy'
    }
    whole_genome_alignment_stats {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/whole_genome_alignments/paftools"
        mode 'copy'
    }
    barrnap {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/barrnap/"
        mode 'copy'
    }
    ssu_16s {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/barrnap/ssu_16s"
        mode 'copy'
    }
    lsu_23s {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/barrnap/lsu_23s"
        mode 'copy'
    }
    lsu_5s {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/barrnap/lsu_5s"
        mode 'copy'
    }
    sina_classified_16s {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/barrnap/ssu_16s/sina/"
        mode 'copy'
    }
    sina_classified_23s {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/barrnap/ssu_23s/sina/"
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
    busco {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/"
        mode 'copy'
    }
    dgenies {
        path "${params.outdir}/${params.sampleid}/${params.taxid}/whole_genome_alignments/dgenies"
        mode 'copy'
    }
}
