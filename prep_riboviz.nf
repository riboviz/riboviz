#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

/*
===================================
RiboViz ribosome profiling workflow
===================================

Running
-------

nextflow run prep_riboviz.nf -params-file <CONFIG>.yaml

where '<CONFIG>' is a YAML configuration file. The YAML configuration
parameters are as follows (all are mandatory unless stated).

Configuration
-------------

Organism data:

* 'orf_fasta_file': Transcript sequences file containing both coding
  regions and flanking regions (FASTA file)
* 'orf_gff_file': Matched genome feature file, specifying coding
  sequences locations (start and stop coordinates) within the
  transcripts (GTF/GFF3 file)
* 'rrna_fasta_file': Ribosomal rRNA and other contaminant sequences to
  avoid aligning to (FASTA file)

Ribosome profiling data:

* 'dir_in': Input directory.
* Either:
  - 'fq_files': Dictionary of FASTQ files to be processed, relative to
     '<dir_in>'. Each item consists of a sample name with a file
     name value (e.g. 'WT3AT: SRR1042864_s1mi.fastq.gz')
* Or:
  - 'multiplex_fq_files': List with a multiplexed FASTQ file,
    relative to '<dir_in>'. If this is provided then the 'fq_files'
    parameter must not be present in the configuration and the
    'sample_sheet' parameter must be present.
  - 'sample_sheet': A sample sheet, relative to '<dir_in>', mandatory
    if 'multiplex_fq_files' is used (tab-separated values file with,
    at least, 'SampleID' and 'TagRead' (barcode) columns)
* If neither or both of 'fq_files' and 'multiplex_fq_files' parameters
  are provided then the workflow will exit.

Indexing:

* 'build_indices': 'TRUE' or 'FALSE', rebuild indices from FASTA
  files? (default 'TRUE')
* 'orf_index_prefix': Prefix for ORF index files, relative to
  '<dir_index>' .
* 'rrna_index_prefix': Prefix for rRNA index files, relative to
  '<dir_index>'.
* 'dir_index': Directory to write indexed files to (default 'index')

Outputs:

* 'dir_tmp': Directory to write temporary files to (default 'tmp')
* 'dir_out': Directory to write temporary files to (default 'output')

Adapter trimming:

* 'adapters': Illumina sequencing adapter(s) to remove.

Barcode and UMI extraction, deduplication, demultiplexing:

* 'extract_umis': 'TRUE' or 'FALSE', extract UMIs after adapter
  trimming? (default 'FALSE')
* 'umi_regexp': UMI-tools-compliant regular expression to extract
  barcodes and UMIs. For details on the regular expression format, see
  UMI-tools documentation on Barcode extraction
  https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction.
  Only required if 'extract_umis' is 'TRUE'.
  - If 'fq_files' are provided then 'umi_regexp' should extract only
    UMIs (i.e. it should contain '<umi>' elements only).
  - If 'multiplex_fq_files' is provided then 'umi_regexp' should
    extract both barcodes and UMIs (i.e. it should contain both
    '<cell>' and '<umi>' elements).
* 'dedup_umis': 'TRUE' or 'FALSE', deduplicate reads using UMI-tools?
  (default 'FALSE')
* 'group_umis': 'TRUE' or 'FALSE', summarise UMI groups both pre- and
  post-deduplication, using UMI-tools? Useful for debugging (default
  'FALSE')
* If 'dedup_umis' is 'TRUE' but 'extract_umis' is 'FALSE' then a
  warning will be displayed, but processing will continue.

Statistics and figure generation input files:

* 'asite_disp_length_file': Summary of read frame displacement from 5'
  end to A-site for each read length based on 'standard' yeast data
  from early ribosome profiling papers (tab-separated values file with
  'read_length', 'asite_disp' columns)
* 'codon_positions_file': Position of codons within each gene (RData
  file)
* 'features_file': Features to correlate with ORFs (tab-separated
  values file with 'ORF', 'Length_log10', 'uATGs', 'FE_atg', 'FE_cap',
  'utr', 'utr_gc', 'polyA' columns)
* 't_rna_file': tRNA estimates file (tab-separated values file with
  'AA', 'Codon', 'tRNA', 'tAI', 'Microarray', 'RNA.seq' columns)

Statistics and figure generation parameters:

* 'buffer': Length of flanking region around the CDS (default 250)
* 'count_reads': 'TRUE' or 'FALSE', scan input, temporary and output
  files and produce counts of reads in each FASTQ, SAM, and BAM file
  processed? (default: 'TRUE')
* 'count_threshold': Remove genes with a read count below this
  threshold, when generating statistics and figures (default 1)
* 'dataset': Human-readable name of the dataset (default 'dataset')
* 'do_pos_sp_nt_freq': 'TRUE' or 'FALSE', calculate position-specific
  nucleotide freqeuency? (default 'TRUE')
* 'is_riboviz_gff': 'TRUE' or 'FALSE', does the GFF file contain 3
  elements per gene - UTR5, CDS, and UTR3? (default 'TRUE')
* 'make_bedgraph': 'TRUE' or 'FALSE', output bedgraph data files in
  addition to H5 files? (default 'TRUE')
* 'max_read_length': Maximum read length in H5 output (default 50)
* 'min_read_length': Minimum read length in H5 output (default 10)
* 'primary_id': Primary gene IDs to access the data (YAL001C, YAL003W,
  etc.) (default 'Name')
* 'rpf': 'TRUE' or 'FALSE', is the dataset an RPF or mRNA dataset?
  (default 'TRUE')
* 'secondary_id': Secondary gene IDs to access the data (COX1, EFB1,
   etc. or 'NULL') (default 'NULL')
* 'stop_in_cds': 'TRUE' or 'FALSE', are stop codons part of the CDS
  annotations in GFF? (default ('FALSE')

General:

* 'num_processes': Number of processes to parallelize over, used by
  specific steps in the workflow (default 1)

Notes
------

Optional inputs follow pattern
https://github.com/nextflow-io/patterns/blob/master/optional-input.nf.
*/

/*
Initialise and validate configuration.

Initialise optional variables to avoid "WARN: Access to undefined
parameter '<PARAM>'" errors.
*/

params.buffer = 250
params.build_indices = true
params.count_reads = true
params.count_threshold = 1
params.dataset = "dataset"
params.dedup_umis = false
params.dir_index = "index"
params.dir_out = "output"
params.dir_tmp = "tmp"
params.do_pos_sp_nt_freq = true
params.extract_umis = false
params.fq_files = [:]
params.group_umis = false
params.is_riboviz_gff = true
params.make_bedgraph = true
params.max_read_length = 50
params.min_read_length = 10
params.multiplex_fq_files = []
params.num_processes = 1
params.primary_id = "Name"
params.rpf = true
params.secondary_id = "NULL"
params.stop_in_cds = false

if (! params.containsKey('adapters')) {
    exit 1, "Undefined adapters (adapters)"
}
if (! params.containsKey('orf_index_prefix')) {
    exit 1, "Undefined ORF index prefix (orf_index_prefix)"
}
if (! params.containsKey('rrna_index_prefix')) {
    exit 1, "Undefined rRNA index prefix (rrna_index_prefix)"
}
if (params.buffer < 0) {
    exit 1, "CDS flanking region length (buffer) is < 0"
}
if (params.count_threshold < 0) {
    exit 1, "Read count threshold (count_threshold) is < 0"
}
if (params.num_processes < 1) {
    exit 1, "Number of processes (num_processes) is < 1"
}
if (params.min_read_length < 1) {
    exit 1, "Minimum read length in H5 output (min_read_length) is < 1"
}
if (params.max_read_length < 1) {
    exit 1, "Maximum read length in H5 output (max_read_length) is < 1"
}
if (params.max_read_length < params.min_read_length) {
    exit 1, "Maximum read length in H5 output is less than minimum read length (max_read_length)"
}
if (! params.secondary_id) {
    secondary_id = "NULL"
} else {
    secondary_id = params.secondary_id
}
if (params.dedup_umis) {
    if (! params.extract_umis) {
        println("Warning: UMI deduplication was requested (dedup_umi: TRUE) but UMI extraction was not (extract_umis: FALSE)")
    }
}
if (params.extract_umis) {
    if (! params.containsKey('umi_regexp')) {
        exit 1, "Undefined barcode/UMI regular expression (umi_regexp) despite UMI extraction being requested (extract_umis: FALSE)"
    }
}

/*
Validate input files.
*/

num_samples = 0
sample_files = [:]
multiplex_files = [:]
multiplex_sample_sheet = Channel.empty()
is_multiplexed = false
if (! params.containsKey('dir_in')) {
    exit 1, "Input directory (dir_in) is undefined"
}
if ((! params.fq_files) && (! params.multiplex_fq_files)) {
    exit 1, "No sample files (fq_files) or multiplexed files (multiplex_fq_files) are defined"
} else if (params.fq_files && params.multiplex_fq_files) {
    exit 1, "Both sample files (fq_files) and multiplexed files (multiplex_fq_files) are defined"
} else if (params.fq_files) {
    // Filter 'params.fq_files' down to those samples that exist.
    for (entry in params.fq_files) {
        sample_file = file("${params.dir_in}/${entry.value}")
        if (sample_file.exists()) {
            sample_files[entry.key] = sample_file
	    num_samples++
        } else {
            println("No such file ($entry.key): $entry.value")
        }
    }
    if (! sample_files) {
        exit 1, "None of the defined sample files (fq_files) exist"
    }
} else {
    // Filter 'params.multiplex_fq_files' down to those files that exist.
    for (entry in params.multiplex_fq_files) {
        multiplex_file = file("${params.dir_in}/${entry}")
        if (multiplex_file.exists()) {
            // Use file base name as key, ensuring that if file
	    // has extension '.fastq.gz' or '.fq.gz' then both
	    // extensions are removed from the name.
            multiplex_file_name = multiplex_file.baseName
            if (multiplex_file_name.endsWith(".fastq")) {
                multiplex_file_name = multiplex_file_name - '.fastq'
            } else if (multiplex_file_name.endsWith(".fq")) {
                multiplex_file_name = multiplex_file_name - '.fq'
            }
            multiplex_files[multiplex_file_name] = multiplex_file
        } else {
            println("No such file: $entry")
        }
    }
    if (! multiplex_files) {
        exit 1, "None of the defined multiplexed files (multiplex_fq_files) exist"
    }
    if (! params.containsKey('sample_sheet')) {
        exit 1, "Undefined sample sheet (sample_sheet)"
    }
    multiplex_sample_sheet = Channel.fromPath(
        "${params.dir_in}/${params.sample_sheet}",
        checkIfExists: true)
    is_multiplexed = true
}

// Create YAML fragment including 'params.fq_files' and
// 'params.multiplex_fq_files' to serve as a configuration file for
// riboviz.tools.count_reads. There is no way to get the location
// of the YAML configuration file itself (i.e. the value of
// '-param-file' from within Nextflow.
Map data_files_config = [:]
data_files_config.fq_files = params.fq_files
data_files_config.multiplex_fq_files = params.multiplex_fq_files
data_files_config_yaml = new Yaml().dump(data_files_config)

// Non-sample-specific input files.
if (! params.containsKey('rrna_fasta_file')) {
    exit 1, "Undefined rRNA FASTA file (rrna_fasta_file)"
}
rrna_fasta = Channel.fromPath(params.rrna_fasta_file,
                              checkIfExists: true)
if (! params.containsKey('orf_fasta_file')) {
    exit 1, "Undefined ORF FASTA file (orf_fasta_file)"
}
orf_fasta = Channel.fromPath(params.orf_fasta_file,
                             checkIfExists: true)
if (! params.containsKey('orf_gff_file')) {
    exit 1, "Undefined ORF GFF file (orf_gff_file)"
}
orf_gff = Channel.fromPath(params.orf_gff_file,
                           checkIfExists: true)

// Optional inputs for generate_stats_figs.R.
// If an optional file is not provided then a 'Missing_<PARAM>' file
// (for example 'Missing_features_file') is created within the 'work/'
// directories for the generateStatsFigs process. This symbolically
// links to a non-existent 'Missing_<PARAM>' file in the users current
// directory. This is not an issue since the files will not be passed
// onto generate_stats_figs.R and no attempt is made to use them. They
// are a side-effect of using the Nextflow pattern for optional
// inputs.
if (params.containsKey('t_rna_file')
    && params.containsKey('codon_positions_file')) {
    t_rna_file = Channel.fromPath(params.t_rna_file,
                                  checkIfExists: true)
    codon_positions_file = Channel.fromPath(params.codon_positions_file, 
                                            checkIfExists: true)
    is_t_rna_and_codon_positions_file = true
} else if ((! params.containsKey('t_rna_file'))
           && (! params.containsKey('codon_positions_file'))) {
    t_rna_file = file("Missing_t_rna_file")
    codon_positions_file = file("Missing_codon_positions_file")
    is_t_rna_and_codon_positions_file = false
} else {
    exit 1, "Either both tRNA estimates (t_rna_file) and codon positions (codon_positions_file) must be defined or neither must be defined"
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
Indexing.
*/ 

// Split channels for use in multiple downstream processes.
orf_fasta.into { orf_fasta_index; orf_fasta_generate_stats_figs }
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

/*
Sample file (fq_files)-specific processes.
*/

process cutAdapters {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(sample_file) \
            from sample_files.collect{ id, file -> [id, file] }
    output:
        tuple val(sample_id), file("trim.fq") into cut_samples
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o trim.fq ${sample_file} -j 0
        """
}

// Route 'cut_samples' channel outputs depending on whether UMIs are
// to be extracted or not.
cut_samples.branch {
    umi_samples: params.extract_umis
    non_umi_samples: ! params.extract_umis
}
.set { cut_samples_branch }

process extractUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(sample_file) \
            from cut_samples_branch.umi_samples
    output:
        tuple val(sample_id), file("extract_trim.fq") \
            into umi_extracted_samples
    when:
        params.extract_umis
    shell:
        """
        umi_tools extract -I ${sample_file} \
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
    publishDir "${params.dir_tmp}", mode: 'copy', overwrite: true
    input:
        tuple val(multiplex_id), file(multiplex_file) \
            from multiplex_files.collect{ id, file -> [id, file] }
    output:
        tuple val(multiplex_id), file("${multiplex_id}_trim.fq") \
            into cut_multiplex
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o ${multiplex_id}_trim.fq ${multiplex_file} -j 0
        """
}

// Route 'cut_multiplex' channel outputs depending on whether UMIs are
// to be extracted or not.
cut_multiplex.branch {
    umi_multiplex: params.extract_umis
    non_umi_multiplex: ! params.extract_umis
}
.set { cut_multiplex_branch }

process extractUmisMultiplex {
    tag "${multiplex_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}", mode: 'copy', overwrite: true
    input:
        tuple val(multiplex_id), file(multiplex_file) \
            from cut_multiplex_branch.umi_multiplex
    output:
        tuple val(multiplex_id), file("${multiplex_id}_extract_trim.fq") \
            into umi_extracted_multiplex
    when:
        params.extract_umis
    shell:
        """
        umi_tools extract -I ${multiplex_file} \
            --bc-pattern="${params.umi_regexp}" \
            --extract-method=regex -S ${multiplex_id}_extract_trim.fq
        """
}

// Combine channels for downstream processing. By definition of
// 'cut_multiplex.branch' only one of the input channels will have
// content.
trimmed_multiplex = cut_multiplex_branch.non_umi_multiplex.mix(
    umi_extracted_multiplex)

// Split channel for use in multiple downstream processes.
multiplex_sample_sheet.into {
    multiplex_sample_sheet_report; multiplex_sample_sheet_demultiplex
}

process demultiplex {
    tag "${multiplex_id}"
    publishDir "${params.dir_tmp}/${multiplex_id}_deplex", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        // Use '.toString' to prevent changing hashes of
	// 'workflow.projectDir' triggering reexecution of this
	// process if 'nextflow run' is run with '-resume'.
        env PYTHONPATH from workflow.projectDir.toString()
        tuple val(multiplex_id), file(multiplex_file) from trimmed_multiplex
        each file(sample_sheet) from multiplex_sample_sheet_demultiplex
    output:
        tuple val(multiplex_id), file("num_reads.tsv") \
                into demultiplexed_num_reads_tsv
        file("*.f*") into demultiplexed_output_fq
    shell:
        """
        python -m riboviz.tools.demultiplex_fastq \
            -1 ${multiplex_file} -s ${sample_sheet} -o . -m 2
        """
}

// 'demultiplexed_output_fq' outputs a single list with all the output
// files. Extract sample IDs from file basenames, filter out
// 'Unassigned' and output tuples of sample IDs and file names as
// separate items onto a new channel.
demultiplexed_output_fq
    .flatten()
    // Use file basename as sample ID.
    .map { [it.baseName, it] }
    // If file was '.fastq|fq.gz' then basename will include
    // '.fq|fastq' so strip that off too.
    .map { n, f -> [n.endsWith(".fq") ? n - ".fq" : n, f] }
    .map { n, f -> [n.endsWith(".fastq") ? n - ".fastq" : n, f] }
    .filter { n, f -> n != "Unassigned" }
    .into { demultiplexed_samples_report; demultiplexed_samples_fq }

if (is_multiplexed) {

    demultiplexed_sample_ids = demultiplexed_samples_report
        .map { n, f -> n }
        .toList() // [] if none
        .view { "Demultiplexed samples: ${it}"}
        // Wrap list in list, so 'merge' below doesn't append lists
        .map { it -> [it] }

    multiplexed_sample_ids = multiplex_sample_sheet_report
        // Extract original sample IDs from sample sheet
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.SampleID } // No output if no 'SampleID' column
        .toList() // [] if no 'SampleID' column
        // Wrap list in list, so 'merge' below doesn't append lists
        .map { it -> [it] }

    multiplexed_sample_ids
        .merge(demultiplexed_sample_ids)
        .map { a, b -> a - b}
        .view { "Non-demultiplexed samples: ${it}" }
}

/*
Sample-specific processes.

Common to both sample files (fq_files) and demultiplexed files.
*/

// Combine channels for downstream processing. By definition of
// 'cut_samples.branch' only one of the input channels will have
// content.
trimmed_samples = cut_samples_branch.non_umi_samples.mix(
    umi_extracted_samples).mix(demultiplexed_samples_fq)

process hisat2rRNA {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
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
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
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
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        // Use '.toString' to prevent changing hashes of
	// 'workflow.projectDir' triggering reexecution of this
	// process if 'nextflow run' is run with '-resume'.
        env PYTHONPATH from workflow.projectDir.toString()
        tuple val(sample_id), file(sam) from orf_map_sams
    output:
        tuple val(sample_id), file("orf_map_clean.sam") \
            into clean_orf_map_sams
        tuple val(sample_id), file("trim_5p_mismatch.tsv") \
            into trim_summary_tsvs
    shell:
        """
        python -m riboviz.tools.trim_5p_mismatch -m 2 \
            -i ${sam} -o orf_map_clean.sam -s trim_5p_mismatch.tsv
        """
}

process samViewSort {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sam) from clean_orf_map_sams
    output:
        tuple val(sample_id), file("orf_map_clean.bam"), \
            file("orf_map_clean.bam.bai") into orf_map_bams
    shell:
        """
        samtools --version
        samtools view -b ${sam} | samtools sort \
            -@ ${params.num_processes} -O bam -o orf_map_clean.bam -
        samtools index orf_map_clean.bam
        """
}

// Route "orf_map_bams" channel outputs depending on whether UMIs are
// to be deduplicated or not.
orf_map_bams.branch {
    umi_bams: params.dedup_umis
    umi_free_bams: ! params.dedup_umis
}
.set { orf_map_bams_branch }

// Split channel for use in multiple downstream processes.
orf_map_bams_branch.umi_bams.into {
    pre_dedup_group_bams; pre_dedup_bams
}

process groupUmisPreDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(bam), file(bam_bai) \
            from pre_dedup_group_bams
    output:
        tuple val(sample_id), file("pre_dedup_groups.tsv") \
            into pre_dedup_groups_tsv
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
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(bam), file(bam_bai) \
            from pre_dedup_bams
    output:
        tuple val(sample_id), file("dedup.bam"), \
            file("dedup.bam.bai") into dedup_bams
        tuple val(sample_id), file("dedup_stats*.tsv") \
            into dedup_stats_tsvs
    when:
        params.dedup_umis
    shell:
        """
        umi_tools dedup -I ${bam} -S dedup.bam --output-stats=dedup_stats
        samtools --version
        samtools index dedup.bam
        """
}

// Split channel for use in multiple downstream processes.
dedup_bams.into { post_dedup_group_bams; post_dedup_bams }

process groupUmisPostDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(bam), file(bam_bai) \
            from post_dedup_group_bams
    output:
        tuple val(sample_id), file("post_dedup_groups.tsv") \
            into post_dedup_groups_tsv
    when:
        params.dedup_umis && params.group_umis
    shell:
        """
        umi_tools group -I ${bam} --group-out post_dedup_groups.tsv
        """
}

// Combine channels for downstream processing. By definition of
// 'orf_map_bams_branch' only one of the input channels will have
// content.
pre_output_bams = orf_map_bams_branch.umi_free_bams.mix(post_dedup_bams)

process outputBams {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(bam), file(bam_bai) \
            from pre_output_bams
    output:
        tuple val(sample_id), file("${sample_id}.bam"), \
            file("${sample_id}.bam.bai") into output_bams
    shell:
        """
        cp ${bam} ${sample_id}.bam
        cp ${bam_bai} ${sample_id}.bam.bai
        """
}

// Split channel for use in multiple downstream processes.
output_bams.into { bedgraph_bams; summary_bams }

process makeBedgraphs {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(bam), file(bam_bai) \
            from bedgraph_bams
    output:
        tuple val(sample_id), file("plus.bedgraph"), \
            file("minus.bedgraph") into bedgraphs
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
    publishDir "${params.dir_out}/${sample_id}", \
        mode: 'copy', overwrite: true
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
    publishDir "${params.dir_out}/${sample_id}", \
        mode: 'copy', overwrite: true
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
        val sample_id into processed_samples
        tuple val(sample_id), file("tpms.tsv") into tpms_tsv
        tuple val(sample_id), file("3nt_periodicity.pdf") \
            into nt3_periodicity_pdf
        tuple val(sample_id), file("3nt_periodicity.tsv") \
            into nt3_periodicity_tsv
        tuple val(sample_id), file("pos_sp_nt_freq.tsv") \
            into pos_sp_nt_freq_tsv
        tuple val(sample_id), file("pos_sp_rpf_norm_reads.pdf") \
            into pos_sp_rpf_norm_reads_pdf
        tuple val(sample_id), file("pos_sp_rpf_norm_reads.tsv") \
            into pos_sp_rpf_norm_reads_tsv
        tuple val(sample_id), file("read_lengths.pdf") \
            into read_lengths_pdf
        tuple val(sample_id), file("read_lengths.tsv") \
            into read_lengths_tsv
        tuple val(sample_id), file("startcodon_ribogridbar.pdf") \
            into start_codon_ribogridbar_pdf
        tuple val(sample_id), file("startcodon_ribogrid.pdf") \
            into start_codon_ribogrid_pdf
        tuple val(sample_id), file("codon_ribodens.pdf") \
            optional (! is_t_rna_and_codon_positions_file) \
            into codon_ribodens_pdf
        tuple val(sample_id), file("codon_ribodens.tsv") \
            optional (! is_t_rna_and_codon_positions_file) \
            into codon_ribodens_tsv
        tuple val(sample_id), file("features.pdf") \
            optional (! is_features_file) into features_tsv
        tuple val(sample_id), file("3ntframe_bygene.tsv") \
            optional (! is_asite_disp_length_file) \
            into nt3frame_bygene_tsv
        tuple val(sample_id), file("3ntframe_propbygene.pdf") \
            optional (! is_asite_disp_length_file) \
            into nt3frame_propbygene_pdf
    shell:
        t_rna_flag = is_t_rna_and_codon_positions_file \
            ? "--t-rna-file=${t_rna}" : ''
        codon_positions_flag = is_t_rna_and_codon_positions_file \
            ? "--codon-positions-file=${codon_positions}" : ''
        features_flag = is_features_file \
            ? "--features-file=${features}" : ''
        asite_disp_length_flag = is_asite_disp_length_file \
            ? "--asite-disp-length-file=${asite_disp_length}" : ''
        count_threshold_flag = params.containsKey('count_threshold') \
            ? "--count-threshold=${params['count_threshold']}": ''
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

processed_samples.view { "Finished processing sample: ${it}" }

// Prefix sample-specific TPMs files, tpms.tsv, with sample ID so all 
// sample-specific TPMs files can be staged into the same directory
// for running collateTpms.
process renameTpms {
    tag "${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(tsv) from tpms_tsv
    output:
        val(sample_id) into sample_tpms_id
        file "${sample_id}_tpms.tsv" into sample_tpms_tsv
    shell:
        """
        cp ${tsv} ${sample_id}_tpms.tsv
        """
}

process collateTpms {
    tag "${sample_ids.join(', ')}"
    publishDir "${params.dir_out}", mode: 'copy', overwrite: true
    input:
        val sample_ids from sample_tpms_id.collect()
        file tsvs from sample_tpms_tsv.collect()
    output:
        file "TPMs_collated.tsv" into collated_tpms_tsv
        val sample_ids into collated_tpms_samples
    shell:
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/collate_tpms.R \
            --sample-subdirs=False \
            --output-dir=. \
            --tpms-file=TPMs_collated.tsv \
            ${sample_ids.join(' ')}
        """
}

process countReads {
    publishDir "${params.dir_out}", mode: 'copy', overwrite: true
    input:
        // Use '.toString' to prevent changing hashes of
	// 'workflow.projectDir' triggering reexecution of this
	// process if 'nextflow run' is run with '-resume'.
        env PYTHONPATH from workflow.projectDir.toString()
        val data_files_config_yaml from data_files_config_yaml
        // Force dependency on output of 'collateTpms' so this process
        // is only run when all other processing has completed.
        val collated_tpms_samples from collated_tpms_samples
    output:
        file "read_counts.tsv" into read_counts_tsv
    when:
        params.count_reads
    shell:
        // 'workflow.projectDir' is directory into which outputs
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
        echo "${data_files_config_yaml}" > data_files.yaml
        python -m riboviz.tools.count_reads \
           -c data_files.yaml \
           -i ${workflow.projectDir}/${params.dir_in} \
           -t ${workflow.projectDir}/${params.dir_tmp} \
           -o ${workflow.projectDir}/${params.dir_out} \
           -r read_counts.tsv  
        """
}

read_counts_tsv.view { "Finished!" }
