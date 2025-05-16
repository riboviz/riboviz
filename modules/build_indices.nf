#!/usr/bin/env nextflow

/*
Indexing.
*/

process buildIndicesrRNA {
    tag "${params.rrna_index_prefix}"
    publishDir "${params.dir_index_env}", mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        path rrna_fasta
    output:
        path "${params.rrna_index_prefix}.*.ht2", emit: built_rrna_index_ht2
    shell:
        """
        hisat2-build --version
        hisat2-build ${rrna_fasta} ${params.rrna_index_prefix}
        """
}

process buildIndicesORF {
    tag "${params.orf_index_prefix}"
    publishDir "${params.dir_index_env}", mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        path orf_fasta
    output:
        path "${params.orf_index_prefix}.*.ht2", emit: built_orf_index_ht2
    shell:
        """
        hisat2-build --version
        hisat2-build ${orf_fasta} ${params.orf_index_prefix}
        """
}