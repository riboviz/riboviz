#!/usr/bin/env nextflow

process trim5pMismatches {
    tag "${sample_id}"
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    errorStrategy 'ignore'
    input:
        env PYTHONPATH
        tuple val(sample_id), path(sample_sam)
    output:
        tuple val(sample_id), path("orf_map_clean.sam"), emit: trim_orf_map_sam
        tuple val(sample_id), path("trim_5p_mismatch.tsv"), emit: trim_summary_tsv
    shell:
        """
        python -m riboviz.tools.trim_5p_mismatch -m 2 \
            -i ${sample_sam} -o orf_map_clean.sam -s trim_5p_mismatch.tsv
        """
}

process samViewSort {
    tag "${sample_id}"
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), path(sample_sam)
    output:
        tuple val(sample_id), path("orf_map_clean.bam"), path("orf_map_clean.bam.bai"), emit: orf_map_bam
    shell:
        memory = params.samsort_memory != null ? "-m ${params.samsort_memory}" : ""
        """
        samtools --version
        samtools view -b ${sample_sam} | samtools sort ${memory} \
            -@ ${params.num_processes} -O bam -o orf_map_clean.bam -
        samtools index orf_map_clean.bam
        """
}



process groupUmisPreDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        tuple val(sample_id), path(sample_bam), path(sample_bam_bai)
    output:
        tuple val(sample_id), path("pre_dedup_groups.tsv") \
            , emit: pre_dedup_group_tsv
    shell:
        """
        umi_tools group -I ${sample_bam} --group-out pre_dedup_groups.tsv
        """
}

process dedupUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        tuple val(sample_id), path(sample_bam), path(sample_bam_bai) 
    output:
        tuple val(sample_id), path("dedup.bam"), \
            path("dedup.bam.bai"), emit: dedup_bam
        tuple val(sample_id), path("dedup_stats*.tsv") \
            , emit: dedup_stats_tsv, optional: (! params.dedup_stats) 
    shell:
        output_stats_flag = params.dedup_stats \
            ? "--output-stats=dedup_stats" : ''
        """
        umi_tools dedup -I ${sample_bam} -S dedup.bam ${output_stats_flag}
        samtools --version
        samtools index dedup.bam
        """
}



process groupUmisPostDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp_env}/${sample_id}", \
        mode: "${params.publish_index_tmp_type}", overwrite: true
    input:
        tuple val(sample_id), path(sample_bam), path(sample_bam_bai)
    output:
        tuple val(sample_id), path("post_dedup_groups.tsv"), emit: post_dedup_group_tsv
    shell:
        """
        umi_tools group -I ${sample_bam} --group-out post_dedup_groups.tsv
        """
}


process outputBams {
    tag "${sample_id}"
    publishDir "${params.dir_out_env}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), path(sample_bam), path(sample_bam_bai)
    output:
        tuple val(sample_id), path("${sample_id}.bam"), \
            path("${sample_id}.bam.bai"), emit: output_bam
    shell:
        """
        cp ${sample_bam} ${sample_id}.bam
        cp ${sample_bam_bai} ${sample_id}.bam.bai
        """
}



process makeBedgraphs {
    tag "${sample_id}"
    publishDir "${params.dir_out_env}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), path(sample_bam), path(sample_bam_bai)
    output:
        tuple val(sample_id), path("plus.bedgraph"), \
            path("minus.bedgraph"), emit: bedgraph
    shell:
        """
        bedtools --version
        bedtools genomecov -ibam ${sample_bam} -trackline -bga -5 \
            -strand + > plus.bedgraph
        bedtools genomecov -ibam ${sample_bam} -trackline -bga -5 \
            -strand - > minus.bedgraph
        """
}

process bamToH5 {
    tag "${sample_id}"
    publishDir "${params.dir_out_env}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), path(sample_bam), \
            path(sample_bam_bai)
        each path(orf_gff)
    output:
        tuple val(sample_id), path("${sample_id}.h5"), path("${sample_id}.h5.*"), emit: h5s
    shell:
        secondary_id_flag = (params.secondary_id != null) \
            ? "--secondary-id=${params.secondary_id}" : ''
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/bam_to_h5.R \
           --num-processes=${params.num_processes} \
           --min-read-length=${params.min_read_length} \
           --max-read-length=${params.max_read_length} \
           --buffer=${params.buffer} \
           --primary-id=${params.primary_id} \
           ${secondary_id_flag} \
           --dataset=${params.dataset} \
           --bam-file=${sample_bam} \
           --hd-file=${sample_id}.h5 \
           --orf-gff-file=${orf_gff} \
           --is-riboviz-gff=${params.is_riboviz_gff} \
           --feature=${params.feature} \
           --stop-in-feature=${params.stop_in_feature}
        """
}