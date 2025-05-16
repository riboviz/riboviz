#!/usr/bin/env nextflow


process hisat2rRNA {
    tag "${sample_id}"
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), path(sample_fq)
        each path(rrna_index_ht2)
    output:
        tuple val(sample_id), path("nonrRNA.fq"), emit: non_rrna_fq
        tuple val(sample_id), path("rRNA_map.sam"), emit: rrna_map_sam
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -N 1 -k 1 \
            --un nonrRNA.fq --no-unal \
            -x ${params.rrna_index_prefix} \
            -S rRNA_map.sam -U ${sample_fq}
        """
}

process hisat2ORF {
    tag "${sample_id}"
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), path(sample_fq)
        each path(orf_index_ht2)
    output:
        tuple val(sample_id), path("unaligned.fq"), emit: unaligned_fq
        tuple val(sample_id), path("orf_map.sam"), emit: trim_5p_mismatches
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} ${params.hisat2_orf_params} \
            --un unaligned.fq -x ${params.orf_index_prefix} \
            -S orf_map.sam -U ${sample_fq}
        """
}