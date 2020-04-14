#!/usr/bin/env nextflow

// Initialise optional variables to avoid "WARN: Access to undefined
// parameter `<PARAM>`" when running Nextflow.
params.make_bedgraph = true
params.extract_umis = false
params.dedup_umis = false
params.group_umis = false

rrna_fasta = Channel.fromPath(params.rrna_fasta_file,
                              checkIfExists: true)

orf_fasta = Channel.fromPath(params.orf_fasta_file,
                             checkIfExists: true)

// Create list of samples whose files exist.
samples = []
for (entry in params.fq_files) {
    sample_file = file("${params.dir_in}/${entry.value}")
    if (sample_file.exists()) {
        samples.add([entry.key, sample_file])
    } else {
        println "Missing file ($entry.key): $entry.value"
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
        file fasta from orf_fasta
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
        tuple val(sample_id), file("trim.fq") into cuts
    shell:
        // TODO configure -j 0 in a more Nextflow-esque way.
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} -o trim.fq ${sample_file} -j 0
        """
}

// Route "cuts" output channel depending on whether UMIs are to be
// extracted.
cuts.branch {
    adapters_umis: params.extract_umis
    adapters: ! params.extract_umis
}
.set { cuts_branch }

process extractUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(sample_file) from cuts_branch.adapters_umis
    output:
        tuple val(sample_id), file("extract_trim.fq") into extract_adapters
    when:
        params.extract_umis
    shell:
        """
        umi_tools extract -I ${sample_file} --bc-pattern="${params.umi_regexp}" --extract-method=regex -S extract_trim.fq
        """
}

// Combine "cuts_branch.adapters" and "extract_adapters" channels and
// route to downstream processing steps. By definition of
// "cuts.branch" only one of the input channels will have content.
trimmed = cuts_branch.adapters.mix(extract_adapters)

process hisat2rRNA {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(fastq) from trimmed
        each file(indices) from rrna_indices
    output:
        tuple val(sample_id), file("nonrRNA.fq") into non_rrnas
        tuple val(sample_id), file("rRNA_map.sam") into rrna_maps
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
        tuple val(sample_id), file(fastq) from non_rrnas
        each file(indices) from orf_indices
    output:
        tuple val(sample_id), file("unaligned.fq") into unaligneds
        tuple val(sample_id), file("orf_map.sam") into orf_maps
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
        tuple val(sample_id), file(sam) from orf_maps
    output:
        tuple val(sample_id), file("orf_map_clean.sam") into clean_orf_maps
        tuple val(sample_id), file("trim_5p_mismatch.tsv") into trim_summaries
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
        tuple val(sample_id), file(sam) from clean_orf_maps
    output:
        tuple val(sample_id), file("orf_map_clean.bam"), file("orf_map_clean.bam.bai") into bams
    shell:
        """
        samtools --version
        samtools view -b ${sam} | samtools sort -@ ${params.num_processes} -O bam -o orf_map_clean.bam -
        samtools index orf_map_clean.bam
        """
}

// Route "bams" output channel depending on whether UMIs are to be
// deduplicated.
bams.branch {
    umi_bams: params.dedup_umis
    umi_free_bams: ! params.dedup_umis
}
.set { bams_branch }

process dedupUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from bams_branch.umi_bams
    output:
        tuple val(sample_id), file("dedup.bam"), file("dedup.bam.bai") into dedup_bams
        tuple val(sample_id), file("dedup_stats*.tsv") into dedup_stats_tsv
    when:
        params.dedup_umis
    shell:
        """
        umi_tools dedup -I ${bam} -S dedup.bam --output-stats=dedup_stats
        samtools --version
        samtools index dedup.bam
        """
}

// Combine "bams_branch.umi_free_bams" and "dedup_bams" channels
// and route to downstream processing steps. By definition of
// "bams_branch" only one of the input channels will have content.
pre_output_bams = bams_branch.umi_free_bams.mix(dedup_bams)

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
output_bams.into { bams_bedgraphs; bams_summary }

process makeBedgraphs {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(bam), file(bam_bai) from bams_bedgraphs
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

// Collect sample IDs and print list.
// Could be used as a model for implementing collect_tpms.R task.
process summarise {
    input:
        val(samples) from bams_summary.map({name, f1, f2 -> return (name) }).collect()
    output:
        val(samples) into summary
    shell:
        """
        echo "Processed samples: ${samples.join(' ')}"
        """
}

// Print output from summarise.
summary.subscribe { println "Processed samples: ${it.join(' ')} "}
