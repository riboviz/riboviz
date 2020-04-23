#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

/*
 * Optional inputs follow pattern
 * https://github.com/nextflow-io/patterns/blob/master/optional-input.nf
 */

/*
 * Validate configuration.
 */
// Initialise optional variables to avoid "WARN: Access to undefined
// parameter `<PARAM>`" when running Nextflow.
params.fq_files = [:]
params.multiplex_fq_files = []
params.make_bedgraph = true
params.extract_umis = false
params.dedup_umis = false
params.group_umis = false
params.count_reads = true

if (! params.secondary_id) {
    secondary_id = "NULL"
} else {
    secondary_id = params.secondary_id
}

if (params.dedup_umis) {
    if (! params.extract_umis) {
        println("WARNING: dedup_umis was TRUE but extract_umis was FALSE")
    }
}

// Filter params.fq_files down to those samples that exist and
// create paths to these files.
sample_files = []
missing = []
for (entry in params.fq_files) {
    sample_file = file("${params.dir_in}/${entry.value}")
    if (sample_file.exists()) {
        sample_files.add([entry.key, sample_file])
    } else {
        println("WARNING: Missing file ($entry.key): $entry.value")
	missing.add(entry)
    }
}
params.fq_files.removeAll{ it -> missing.contains(it) }

// Filter params.multiplex_fq_files down to those files that exist
// andcreate paths to these files.

multiplex_sample_files = []
missing = []
for (entry in params.multiplex_fq_files) {
    sample_file = file("${params.dir_in}/${entry}")
    if (sample_file.exists()) {
        multiplex_sample_files.add(sample_file)
    } else {
        println("WARNING: Missing multiplexed file: $entry")
	missing.add(entry)
    }
}
params.multiplex_fq_files.removeAll{ it -> missing.contains(it) }

/*
 * Set up non-sample-specific input files.
 */
rrna_fasta = Channel.fromPath(params.rrna_fasta_file,
                              checkIfExists: true)
orf_gff = Channel.fromPath(params.orf_gff_file,
                           checkIfExists: true)
orf_fasta = Channel.fromPath(params.orf_fasta_file,
                             checkIfExists: true)
/*
 * Set up optional inputs for generate_stats_figs.R.
 *
 * If an optional file is not provided then a "Missing_<PARAM>" file
 * (for example "Missing_features_file") is created within the work/
 * directories for the generateStatsFigs process. This symbolically
 * links to a non-existent "Missing_<PARAM>" file in the users current
 * directory. This is not an issue since the files will not be passed
 * onto generate_stats_figs.R and no attempt is made to use them. They
 * are a side-effect of using the Nextflow pattern for optional inputs,
 * https://github.com/nextflow-io/patterns/blob/master/optional-input.nf.
 */
if (params.containsKey('t_rna_file') && params.containsKey('codon_positions_file')) {
    t_rna_file = Channel.fromPath(params.t_rna_file, checkIfExists: true)
    codon_positions_file = Channel.fromPath(params.codon_positions_file,
                                            checkIfExists: true)
    is_t_rna_and_codon_positions_file = true
} else if ((! params.containsKey('t_rna_file')) && (! params.containsKey('codon_positions_file')))
{
    t_rna_file = file("Missing_t_rna_file")
    codon_positions_file = file("Missing_codon_positions_file")
    is_t_rna_and_codon_positions_file = false
} else {
    error "ERROR: Either both t_rna_file and codon_positions_file must be provided or none must be provided."
}
if (params.containsKey('features_file')) {
    features_file = Channel.fromPath(params.features_file,
                                     checkIfExists: true)
    is_features_file = true
} else {
    features_file = file("Missing_features_file")
    is_features_file = false
}
if (params.containsKey('asite_disp_length_file')) {
    asite_disp_length_file = Channel.fromPath(params.asite_disp_length_file,
                                              checkIfExists: true)
    is_asite_disp_length_file = true
} else {
    asite_disp_length_file = file("Missing_aside_disp_length_file")
    is_asite_disp_length_file = false
}

/*
 * Workflow processes.
 */ 

// Split "orf_fasta" channel so can use as input to multiple
// downstream tasks.
orf_fasta.into { orf_fasta_index; orf_fasta_generate_stats_figs }

// Split "orf_gff" channel so can use as input to multiple
// downstream tasks.
orf_gff.into { orf_gff_bam_to_h5; orf_gff_generate_stats_figs }

process buildIndicesrRNA {
    tag "${params.rrna_index_prefix}"
    publishDir "${params.dir_index}", mode: 'copy', overwrite: true
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
    publishDir "${params.dir_index}", mode: 'copy', overwrite: true
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
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(sample_file) from sample_files
    output:
        tuple val(sample_id), file("trim.fq") into cut_samples
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o trim.fq ${sample_file} -j 0
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
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(sample_file) from cut_samples_branch.umi_samples
    output:
        tuple val(sample_id), file("extract_trim.fq") into umi_extracted_samples
    when:
        params.extract_umis
    shell:
        """
        umi_tools extract -I ${sample_file} \
            --bc-pattern="${params.umi_regexp}" \
            --extract-method=regex -S extract_trim.fq
        """
}

// Combine channels and route to downstream processing steps. By
// definition of "cut_samples.branch" only one of the input channels
// will have content. 
trimmed_samples = cut_samples_branch.non_umi_samples.mix(umi_extracted_samples)

process hisat2rRNA {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
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
        hisat2 -p ${params.num_processes} -N 1 -k 1 \
            --un nonrRNA.fq -x ${params.rrna_index_prefix} \
            -S rRNA_map.sam -U ${fastq}
        """
}

process hisat2ORF {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
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
        hisat2 -p ${params.num_processes} -k 2 \
            --no-spliced-alignment --rna-strandness F --no-unal \
            --un unaligned.fq -x ${params.orf_index_prefix} \
            -S orf_map.sam -U ${fastq}
        """
}

process trim5pMismatches {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        env PYTHONPATH from workflow.projectDir
        tuple val(sample_id), file(sam) from orf_map_sams
    output:
        tuple val(sample_id), file("orf_map_clean.sam") into clean_orf_map_sams
        tuple val(sample_id), file("trim_5p_mismatch.tsv") into trim_summary_tsvs
    shell:
        """
        python -m riboviz.tools.trim_5p_mismatch -m 2 \
            -i ${sam} -o orf_map_clean.sam -s trim_5p_mismatch.tsv
        """
}

process samViewSort {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sam) from clean_orf_map_sams
    output:
        tuple val(sample_id), file("orf_map_clean.bam"), file("orf_map_clean.bam.bai") into orf_map_bams
    shell:
        """
        samtools --version
        samtools view -b ${sam} | samtools sort \
            -@ ${params.num_processes} -O bam -o orf_map_clean.bam -
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
orf_map_bams_branch.umi_bams.into {
    pre_dedup_group_bams; pre_dedup_bams
}

process groupUmisPreDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from pre_dedup_group_bams
    output:
        tuple val(sample_id), file("pre_dedup_groups.tsv") into pre_dedup_groups_tsv
    when:
        params.dedup_umis && params.group_umis
    shell:
        """
        umi_tools group -I ${bam} --group-out pre_dedup_groups.tsv
        """
}

process dedupUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
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
    publishDir "${params.dir_tmp}/${sample_id}", mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from post_dedup_group_bams
    output:
        tuple val(sample_id), file("post_dedup_groups.tsv") into post_dedup_groups_tsv
    when:
        params.dedup_umis && params.group_umis
    shell:
        """
        umi_tools group -I ${bam} --group-out post_dedup_groups.tsv
        """
}

// Combine "orf_map_bams_branch.umi_free_bams" and "post_dedup_bams"
// channels and route to downstream processing steps. By definition of
// "orf_map_bams_branch" only one of the input channels will have
// content.
pre_output_bams = orf_map_bams_branch.umi_free_bams.mix(post_dedup_bams)

process outputBams {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", mode: 'copy', overwrite: true
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
    publishDir "${params.dir_out}/${sample_id}", mode: 'copy', overwrite: true
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
        bedtools genomecov -ibam ${bam} -trackline -bga -5 \
            -strand + > plus.bedgraph
        bedtools genomecov -ibam ${bam} -trackline -bga -5 \
            -strand - > minus.bedgraph
        """
}

process bamToH5 {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from summary_bams
        each file(gff) from orf_gff_bam_to_h5
    output:
        tuple val(sample_id), file("${sample_id}.h5") into h5s
    shell:
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/bam_to_h5.R \
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

process generateStatsFigs {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(h5) from h5s
	each file(fasta) from orf_fasta_generate_stats_figs
        each file(gff) from orf_gff_generate_stats_figs
	each file(t_rna) from t_rna_file
	each file(codon_positions) from codon_positions_file
	each file(features) from features_file
	each file(asite_disp_length) from asite_disp_length_file
    output:
        tuple val(sample_id), file("tpms.tsv") into tpms_tsv
        tuple val(sample_id), file("*.pdf") into stats_figs_pdfs
        tuple val(sample_id), file("*.tsv") into stats_figs_tsvs
        tuple val(sample_id), file("codon_ribodens.pdf") optional (! is_t_rna_and_codon_positions_file) into codon_ribodens_pdf
        tuple val(sample_id), file("codon_ribodens.tsv") optional (! is_t_rna_and_codon_positions_file) into codon_ribodens_tsv
        tuple val(sample_id), file("features.pdf") optional (! is_features_file) into features_tsv
        tuple val(sample_id), file("3ntframe_bygene.tsv") optional (! is_asite_disp_length_file) into nt3frame_bygene_tsv
        tuple val(sample_id), file("3ntframe_propbygene.pdf") optional (! is_asite_disp_length_file) into nt3frame_propbygene_pdf
    shell:
        t_rna_flag = is_t_rna_and_codon_positions_file ? "--t-rna-file=${t_rna}" : ''
        codon_positions_flag = is_t_rna_and_codon_positions_file ? "--codon-positions-file=${codon_positions}" : ''
        features_flag = is_features_file ? "--features-file=${features}" : ''
        asite_disp_length_flag = is_asite_disp_length_file ? "--asite-disp-length-file=${asite_disp_length}" : ''
        count_threshold_flag = params.containsKey('count_threshold') ? "--count-threshold=${params['count_threshold']}": ''
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/generate_stats_figs.R \
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
           ${t_rna_flag} \
           ${codon_positions_flag} \
           ${features_flag} \
	   --orf-gff-file=${gff} \
           ${asite_disp_length_flag} \
           ${count_threshold_flag}
        """
}

process prepareCollateTpms {
    tag "${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(tsv) from tpms_tsv
    output:
        val sample_id into sample_tpms_id
        file "${sample_id}_tpms.tsv" into sample_tpms_tsv
    shell:
        """
        cp ${tsv} ${sample_id}_tpms.tsv
        """
}

process collateTpms {
    publishDir "${params.dir_out}", mode: 'copy', overwrite: true
    input:
        val sample_ids from sample_tpms_id.collect()
        file tsvs from sample_tpms_tsv.collect()
    output:
        file "TPMs_collated.tsv" into collated_tpms_tsv
        val sample_ids into completed_samples
    shell:
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/collate_tpms.R \
            --sample-subdirs=False \
            --output-dir=. \
            --tpms-file=TPMs_collated.tsv \
            ${sample_ids.join(' ')}
        """
}

completed_samples.subscribe { println("Processed samples: ${it.join(' ')}") }

// Create YAML fragment including params.fq_files and
// params.multiplex_fq_files to serve as a configuration file for
// riboviz.tools.count_reads. There is no way to get the location
// of the RiboViz YAML configuration file itself from within Nextflow
// (there is no way to access the "-param-file" argument to
// Nextflow).
Map count_reads_config = [:]
count_reads_config.fq_files = params.fq_files
count_reads_config.multiplex_fq_files = params.multiplex_fq_files
count_reads_config_yaml = new Yaml().dump(count_reads_config)

process countReads {
    publishDir "${params.dir_out}", mode: 'copy', overwrite: true
    input:
        env PYTHONPATH from workflow.projectDir
        val count_reads_config_yaml from count_reads_config_yaml
	// Force dependency on output of collateTpms so this process
        // is only run when all other processing has completed.
	val completed_samples from completed_samples
    output:
        file "read_counts.tsv" into read_counts_tsv
    when:
        params.count_reads
    shell:
        // workflow.projectDir is directory into which outputs
	// have been published.
	// TODO: It would be preferable to:
	// 1. Stage these directories into this process's
	//    work directory. If this is possible, and how to do it,
        //    if so, has not been determined.
	// 2. Stage outputs from the processes for which reads
	//    are to be counted into this process's work directory.
        //    This would require a new implementation of
	//    riboviz.tools.count_reads.
        """
        echo "${count_reads_config_yaml}" > sample_files.yaml
        python -m riboviz.tools.count_reads \
           -c sample_files.yaml \
           -i ${workflow.projectDir}/${params.dir_in} \
           -t ${workflow.projectDir}/${params.dir_tmp} \
           -o ${workflow.projectDir}/${params.dir_out} \
           -r read_counts.tsv  
        """
}
