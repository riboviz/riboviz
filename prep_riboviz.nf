#!/usr/bin/env nextflow

// Initialise optional variables to avoid "WARN: Access to undefined
// parameter `<PARAM>`" when running Nextflow.
params.make_bedgraph = true
params.extract_umis = false
params.dedup_umis = false
params.group_umis = false
if (! params.secondary_id) {
    secondary_id = "NULL"
} else {
    secondary_id = params.secondary_id
}
if (params.dedup_umis) {
    if (! params.extract_umis) {
        println("WARNING: dedup_umis was TRUE but extract_umis was FALSE.")
    }
}

rrna_fasta = Channel.fromPath(params.rrna_fasta_file,
                              checkIfExists: true)

orf_fasta = Channel.fromPath(params.orf_fasta_file,
                             checkIfExists: true)
// Split "orf_fasta" channel so can use as input to multiple
// downstream tasks.
orf_fasta.into { orf_fasta_index; orf_fasta_generate_stats_figs }

orf_gff = Channel.fromPath(params.orf_gff_file,
                           checkIfExists: true)
// Split "orf_fasta" channel so can use as input to multiple
// downstream tasks.
orf_gff.into { orf_gff_bam_to_h5; orf_gff_generate_stats_figs }

// Create list of samples whose files exist.
samples = []
for (entry in params.fq_files) {
    sample_file = file("${params.dir_in}/${entry.value}")
    if (sample_file.exists()) {
        samples.add([entry.key, sample_file])
    } else {
        println("WARNING: Missing file ($entry.key): $entry.value")
    }
}

process buildIndicesrRNA {
    tag "${params.rrna_index_prefix}"
    publishDir "${params.dir_index}"
    input:
        file fasta from rrna_fasta
    output:
        file "${params.rrna_index_prefix}.*.ht2" into rrna_indices
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${fasta} ${params.rrna_index_prefix}
        """
}

process buildIndicesORF {
    tag "${params.orf_index_prefix}"
    publishDir "${params.dir_index}"
    input:
        file fasta from orf_fasta_index
    output:
        file "${params.orf_index_prefix}.*.ht2" into orf_indices
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${fasta} ${params.orf_index_prefix}
        """
}

process cutAdapters {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(sample_file) from samples
    output:
        tuple val(sample_id), file("trim.fq") into cut_samples
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} -o trim.fq ${sample_file} -j 0
        """
}

// Route "cut_samples" channel depending on whether UMIs are to be
// extracted or not.
cut_samples.branch {
    umi_samples: params.extract_umis
    non_umi_samples: ! params.extract_umis
}
.set { cut_samples_branch }

process extractUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(sample_file) from cut_samples_branch.umi_samples
    output:
        tuple val(sample_id), file("extract_trim.fq") into umi_extracted_samples
    when:
        params.extract_umis
    shell:
        """
        umi_tools extract -I ${sample_file} --bc-pattern="${params.umi_regexp}" --extract-method=regex -S extract_trim.fq
        """
}

// Combine channels and route to downstream processing steps. By
// definition of "cut_samples.branch" only one of the input channels
// will have content.
trimmed_samples = cut_samples_branch.non_umi_samples.mix(umi_extracted_samples)

process hisat2rRNA {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(fastq) from trimmed_samples
        each file(indices) from rrna_indices
    output:
        tuple val(sample_id), file("nonrRNA.fq") into non_rrna_fqs
        tuple val(sample_id), file("rRNA_map.sam") into rrna_map_sams
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -N 1 -k 1 --un nonrRNA.fq -x ${params.rrna_index_prefix} -S rRNA_map.sam -U ${fastq}
        """
}

process hisat2ORF {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(fastq) from non_rrna_fqs
        each file(indices) from orf_indices
    output:
        tuple val(sample_id), file("unaligned.fq") into unaligned_fqs
        tuple val(sample_id), file("orf_map.sam") into orf_map_sams
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -k 2 --no-spliced-alignment --rna-strandness F --no-unal --un unaligned.fq -x ${params.orf_index_prefix} -S orf_map.sam -U ${fastq}
        """
}

process trim5pMismatches {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sam) from orf_map_sams
    output:
        tuple val(sample_id), file("orf_map_clean.sam") into clean_orf_map_sams
        tuple val(sample_id), file("trim_5p_mismatch.tsv") into trim_summary_tsvs
    shell:
        """
        python -m riboviz.tools.trim_5p_mismatch -m 2 -i ${sam} -o orf_map_clean.sam -s trim_5p_mismatch.tsv
        """
}

process samViewSort {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sam) from clean_orf_map_sams
    output:
        tuple val(sample_id), file("orf_map_clean.bam"), file("orf_map_clean.bam.bai") into orf_map_bams
    shell:
        """
        samtools --version
        samtools view -b ${sam} | samtools sort -@ ${params.num_processes} -O bam -o orf_map_clean.bam -
        samtools index orf_map_clean.bam
        """
}

// Route "orf_map_bams" output channel depending on whether UMIs are
// to be deduplicated or not.
orf_map_bams.branch {
    umi_bams: params.dedup_umis
    umi_free_bams: ! params.dedup_umis
}
.set { orf_map_bams_branch }

// Split "orf_map_bams_branch.umi_bams" channel so can use as input to
// multiple downstream tasks.
orf_map_bams_branch.umi_bams.into { pre_dedup_group_bams; pre_dedup_bams }

process groupUmisPreDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from pre_dedup_group_bams
    output:
        tuple val(sample_id), file("pre_dedup_groups.tsv") into pre_dedup_groups_tsv
    when:
        params.group_umis
    shell:
        """
        umi_tools group -I ${bam} --group-out pre_dedup_groups.tsv
        """
}

process dedupUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from pre_dedup_bams
    output:
        tuple val(sample_id), file("dedup.bam"), file("dedup.bam.bai") into dedup_bams
        tuple val(sample_id), file("dedup_stats*.tsv") into dedup_stats_tsvs
    when:
        params.dedup_umis
    shell:
        """
        umi_tools dedup -I ${bam} -S dedup.bam --output-stats=dedup_stats
        samtools --version
        samtools index dedup.bam
        """
}

// Split "dedup_bams" channel so can use as input to multiple
// downstream tasks.
dedup_bams.into { post_dedup_group_bams; post_dedup_bams }

process groupUmisPostDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from post_dedup_group_bams
    output:
        tuple val(sample_id), file("post_dedup_groups.tsv") into post_dedup_groups_tsv
    when:
        params.group_umis
    shell:
        """
        umi_tools group -I ${bam} --group-out post_dedup_groups.tsv
        """
}

// Combine "orf_map_bams_branch.umi_free_bams" and "post_dedup_bams"
// channels and route to downstream processing steps. By definition of
// "orf_map_bams_branch" only one of the input channels will have content.
pre_output_bams = orf_map_bams_branch.umi_free_bams.mix(post_dedup_bams)

process outputBams {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from pre_output_bams
    output:
        tuple val(sample_id), file("${sample_id}.bam"), file("${sample_id}.bam.bai") into output_bams
    shell:
        """
        cp ${bam} ${sample_id}.bam
        cp ${bam_bai} ${sample_id}.bam.bai
        """
}

// Split "output_bams" channel so can use as input to multiple
// downstream tasks.
output_bams.into { bedgraph_bams; summary_bams }

process makeBedgraphs {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from bedgraph_bams
    output:
        tuple val(sample_id), file("plus.bedgraph"), file("minus.bedgraph") into bedgraphs
    when:
        params.make_bedgraph
    shell:
        """
        bedtools --version
        bedtools genomecov -ibam ${bam} -trackline -bga -5 -strand + > plus.bedgraph
        bedtools genomecov -ibam ${bam} -trackline -bga -5 -strand - > minus.bedgraph
        """
}

process bamToH5 {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from summary_bams
        each file(gff) from orf_gff_bam_to_h5
    output:
        tuple val(sample_id), file("${sample_id}.h5") into h5s
    shell:
        """
        Rscript --vanilla /home/ubuntu/riboviz/rscripts/bam_to_h5.R \
           --num-processes=${params.num_processes} \
           --min-read-length=${params.min_read_length} \
           --max-read-length=${params.max_read_length} \
           --buffer=${params.buffer} \
           --primary-id=${params.primary_id} \
           --secondary-id=${secondary_id} \
           --dataset=${params.dataset} \
           --bam-file=${bam} \
           --hd-file=${sample_id}.h5 \
           --orf-gff-file=${gff} \
           --is-riboviz-gff=${params.is_riboviz_gff} \
           --stop-in-cds=${params.stop_in_cds}
        """
}

t_rna = Channel.fromPath(params.t_rna_file,
                         checkIfExists: true)
codon_positions = Channel.fromPath(params.codon_positions_file,
                                   checkIfExists: true)
features = Channel.fromPath(params.features_file,
                            checkIfExists: true)
asite_disp_length = Channel.fromPath(params.asite_disp_length_file,
                                     checkIfExists: true)

process generateStatsFigs {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(h5) from h5s
	each file(fasta) from orf_fasta_generate_stats_figs
        each file(gff) from orf_gff_generate_stats_figs
	each file(t_rna) from t_rna
	each file(codon_positions) from codon_positions
	each file(features) from features
	each file(asite_disp_length) from asite_disp_length
    output:
        tuple val(sample_id), file("tpms.tsv") into tpms_tsv
        tuple val(sample_id), file("*.pdf") into stats_figs_pdfs
        tuple val(sample_id), file("*.tsv") into stats_figs_tsvs
    shell:
        """
        Rscript --vanilla /home/ubuntu/riboviz/rscripts/generate_stats_figs.R \
           --num-processes=${params.num_processes} \
           --min-read-length=${params.min_read_length} \
           --max-read-length=${params.max_read_length} \
           --buffer=${params.buffer} \
           --primary-id=${params.primary_id} \
           --dataset=${params.dataset} \
           --hd-file=${h5} \
           --orf-fasta-file=${fasta} \
           --rpf=${params.rpf} \
           --output-dir=. \
           --do-pos-sp-nt-freq=${params.do_pos_sp_nt_freq} \
           --t-rna-file=${t_rna} \
           --codon-positions-file=${codon_positions} \
           --features-file=${features} \
	   --orf-gff-file=${gff} \
           --asite-disp-length-file=${asite_disp_length} \
           --count-threshold=${params.count_threshold}
        """
}

process prepareCollateTpms {
    tag "${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(tsv) from tpms_tsv
    output:
        val(sample_id) into sample_tpms_id
        file("${sample_id}_tpms.tsv") into sample_tpms_tsv
    shell:
        """
        cp ${tsv} ${sample_id}_tpms.tsv
        """
}

process collateTpms {
    errorStrategy 'ignore'
    publishDir "${params.dir_out}"
    input:
        val(sample_ids) from sample_tpms_id.collect()
        file(tsvs) from sample_tpms_tsv.collect()
    output:
        file("TPMs_collated.tsv") into collated_tpms_tsv
        val(sample_ids) into completed_samples
    shell:
        """
        Rscript --vanilla /home/ubuntu/riboviz/rscripts/collate_tpms.R \
            --sample-subdirs=False \
            --output-dir=. \
            --tpms-file=TPMs_collated.tsv \
            ${sample_ids.join(' ')}
        """
}

completed_samples.subscribe { println("Processed samples: ${it.join(' ')}") }
