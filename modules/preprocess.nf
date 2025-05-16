#!/usr/bin/env nextflow

/*
Sample file (fq_files)-specific processes.
*/

process cutAdapters {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        tuple val(sample_id), path(sample_fq) 
    output:
        tuple val(sample_id), path("trim.fq"), emit: cut_fq
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o trim.fq ${sample_fq} -j 1
        """
}

process extractUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        tuple val(sample_id), path(sample_fq)
    output:
        tuple val(sample_id), path("extract_trim.fq"), emit: umi_extract_fq
    shell:
        """
        umi_tools extract -I ${sample_fq} \
            --bc-pattern="${params.umi_regexp}" \
            --extract-method=regex -S extract_trim.fq
        """
}

/*
Multiplexed files (multiplex_fq_files)-specific processes.
*/

process cutAdaptersMultiplex {
    tag "${multiplex_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp_env}", mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        tuple val(multiplex_id), path(multiplex_fq)
    output:
        tuple val(multiplex_id), path("${multiplex_id}_trim.fq"), emit: cut_multiplex_fq
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o ${multiplex_id}_trim.fq ${multiplex_fq} -j 0
        """
}


process extractUmisMultiplex {
    tag "${multiplex_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp_env}", mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        tuple val(multiplex_id), path(multiplex_fq)
    output:
        tuple val(multiplex_id), path("${multiplex_id}_extract_trim.fq"), emit: umi_extract_multiplex_fq
    when:
        params.extract_umis && is_multiplexed
    shell:
        """
        umi_tools extract -I ${multiplex_fq} \
            --bc-pattern="${params.umi_regexp}" \
            --extract-method=regex -S ${multiplex_id}_extract_trim.fq
        """
}


process demultiplex {
    tag "${multiplex_id}"
    publishDir "${params.dir_tmp_env}/${multiplex_id}_deplex", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    errorStrategy 'ignore'
    input:
        env PYTHONPATH
        tuple val(multiplex_id), path(multiplex_fq)
        each path(sample_sheet_tsv)
    output:
        tuple val(multiplex_id), path("num_reads.tsv"), emit: demultiplex_num_reads_tsv
        path("*.f*"), emit: demultiplex_fq
    shell:
        """
        python -m riboviz.tools.demultiplex_fastq \
            -1 ${multiplex_fq} -s ${sample_sheet_tsv} -o . -m 2
        """
}