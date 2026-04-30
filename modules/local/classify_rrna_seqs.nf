#!/usr/bin/env nextflow

// Use SINA to align our extracted rRNA sequences to a SILVA SSU (16S) database 
// Note FASTA 16S should be a fasta with a single 16s sequence
process SINA_SEARCH_AND_CLASSIFY {

    container "community.wave.seqera.io/library/sina_trimal:0983aef64094b81d"

    input:
    tuple val(sampleid), val(taxid), path(fasta_16s)
    path silva_arb_database

    output:
    tuple val(sampleid), val(taxid), path("${sampleid}.${taxid}.sina.search.csv"), emit: searchcsv
    tuple val(sampleid), val(taxid), path("${sampleid}.${taxid}.sina.search.fasta"), emit: searchfasta
    tuple val(sampleid), val(taxid), path("${sampleid}.${taxid}.sina.search.trimmed.fasta"), emit: trimmed
    tuple val(sampleid), val(taxid), path("${sampleid}.${taxid}.sina.search.csv"), path("${sampleid}.${taxid}.sina.search.fasta"), path("${sampleid}.${taxid}.sina.search.csv"), emit: all

    script:
    """
    set -euo pipefail
    
    # Ask sina to do 3 operations. 
    # 1. Add the sequence to an existing MSA (the silva arb 16S database)
    # 2. Execute a homology search based on the computed alignment, identifying generating a classification. 
    # 3. Output an fasta MSA with our input sequences + 'neighbour' sequences (based on homology search). 
    # Note the reason that sina doesn't just compute pairwise similarity with every ref and take the minimum is because it would take a very long time

    # Search and classify. Each result sequence must have at least 90% fractional identity with query
    # We will only return 10 search results per query sequence (from the 10 most similar seqs)
    sina --search \
    -i ${fasta_16s} -r ${silva_arb_database} \
    --add-relatives 10 --overhang remove \
    --search-min-sim 0.90 --search-max-result 10 --lca-fields tax_gtdb \
    -o ${sampleid}.${taxid}.sina.search.fasta \
    -o ${sampleid}.${taxid}.sina.search.csv 

    # Since sina adds sequences to a large refrence msa alignment the final MSA can have loads of gaps.
    # We remove columns which represent gaps across all sequences in the msa (input seq + neighbours). 
    # This trimmed msa is what should be input into pairwise-distance calculation and tree-construction tools
    trimAl -noallgaps -in ${sampleid}.${taxid}.sina.search.fasta -out ${sampleid}.${taxid}.sina.search.trimmed.fasta
    """
}
