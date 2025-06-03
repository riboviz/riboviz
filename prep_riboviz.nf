#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml


/*
===================================
riboviz ribosome profiling workflow
===================================
*/

/**
 * Log a help message.
 */
def helpMessage() {
    help = """
    Usage
    -----

    Run:

    \$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE> [--help]

    where:

    * '<CONFIG_FILE>' is a YAML configuration file. The YAML
      configuration parameters are described below (all are mandatory
      unless otherwise stated).
    * '--help' displays this help information and exits.
    * Configuration parameters can also be provided via the
      command-line in the form '--<PARAMETER>=<VALUE>' (for example
      '--make_bedgraph=FALSE').

    To specify values for environment variables (see Environment
    variable and configuration tokens), you have two options, where:

    * '<SAMPLES_DIRECTORY>' is a directory with input files.
    * '<ORGANISMS_DIRECTORY>' is a directory with input files.
    * '<DATA_DIRECTORY>' is a directory with input files.

    The options are:

    1. Specify environment variables with the paths to the directories
       on the same line as your command to run the workflow. The
       values will be used for this run of the workflow only. For
       example:

    \$ RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY> \\
       RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY> \\
       RIBOVIZ_DATA=<DATA_DIRECTORY> \\
       nextflow run prep_riboviz.nf -params-file <CONFIG_FILE>

    2. Define values for the environment variables within your bash
       shell. The values will be available for successive runs of the
       workflow. The values need to be defined using 'export' so they
       are available to 'nextflow' when it runs. For example:

    \$ export RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY>
    \$ export RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY>
    \$ export RIBOVIZ_DATA=<DATA_DIRECTORY>
    \$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE>

    The above approaches can be combined i.e. you can define variables
    using 'export' (2) but provide other values as part of the command
    to run the workflow (1). Values provided within the command take
    precedence over those defined via 'export'.

    Configuration
    -------------

    Organism data:

    * 'orf_fasta_file': Transcript sequences file containing both
      coding regions and flanking regions (FASTA file)
    * 'orf_gff_file': Matched genome feature file, specifying coding
      sequences locations (start and stop coordinates) within the
      transcripts (GTF/GFF3 file)
    * 'rrna_fasta_file': Ribosomal rRNA and other contaminant
      sequences to avoid aligning to (FASTA file)

    Ribosome profiling data:

    * 'dir_in': Input directory.
    * Either:
      - 'fq_files': Dictionary of FASTQ files to be processed,
        relative to '<dir_in>'. Each item consists of a sample name
        with a file name value
        (e.g. 'WT3AT: SRR1042864_s1mi.fastq.gz')
    * Or:
      - 'multiplex_fq_files': List with a multiplexed FASTQ file,
        relative to '<dir_in>'. If this is provided then the
        'fq_files' parameter must not be present in the configuration
        and the 'sample_sheet' parameter must be present.
      - 'sample_sheet': A sample sheet, relative to '<dir_in>',
        mandatory if 'multiplex_fq_files' is used (tab-separated
        values file with, at least, 'SampleID' and 'TagRead' (barcode)
        columns)
    * If neither or both of 'fq_files' and 'multiplex_fq_files'
      parameters are provided then the workflow will exit.

    Indexing:

    * 'build_indices': Rebuild indices from FASTA files? (default
      'TRUE'). If 'FALSE' then 'dir_index' is expected to contain the
      index files.
    * 'orf_index_prefix': Prefix for ORF index files, relative to
      '<dir_index>' .
    * 'rrna_index_prefix': Prefix for rRNA index files, relative to
      '<dir_index>'.
    * 'dir_index': Directory to write indexed files to (default 'index')

    Outputs:

    * 'dir_tmp': Directory to write temporary files to (default 'tmp')
    * 'dir_out': Directory to write temporary files to
      (default 'output')

    Adapter trimming:

    * 'adapters': Illumina sequencing adapter(s) to remove.

    Alignment:

    * 'hisat2_orf_params': Command-line parameters for hisat2
      invocation to align ORFs to index files
      (default "-k 2 --no-spliced-alignment --rna-strandness F --no-unal")

    Barcode and UMI extraction, deduplication, demultiplexing:

    * 'extract_umis': Extract UMIs after adapter trimming? (default
      'FALSE')
    * 'umi_regexp': UMI-tools-compliant regular expression to extract
      barcodes and UMIs. For details on the regular expression format,
      see UMI-tools documentation on Barcode extraction
      https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction.
      Only required if 'extract_umis' is 'TRUE'.
      - If 'fq_files' are provided then 'umi_regexp' should extract
        only UMIs (i.e. it should contain '<umi>' elements only).
      - If 'multiplex_fq_files' is provided then 'umi_regexp' should
        extract both barcodes and UMIs (i.e. it should contain both
        '<cell>' and '<umi>' elements).
    * 'dedup_umis': Deduplicate reads using UMI-tools? (default
      'FALSE')
    * 'dedup_stats': Output UMI deduplication statistics? (default
      'TRUE')
    * 'group_umis': Summarise UMI groups both pre- and
      post-deduplication using UMI-tools? Useful for debugging
      (default 'FALSE')
    * If 'dedup_umis' is 'TRUE' but 'extract_umis' is 'FALSE' then a
      warning will be displayed, but processing will continue.
    * 'trim_5p_mismatches': Trim mismatched 5' base? (default 'TRUE')

    Statistics and figure generation input files:

    * 'asite_disp_length_file': Summary of read frame displacement
      from 5' end to A-site for each read length based on 'standard'
      yeast data from early ribosome profiling papers (tab-separated
      values file with 'read_length', 'asite_disp' columns) (optional)
    * 'codon_positions_file': Position of codons within each gene
      (RData file) (optional)
    * 'features_file': Features to correlate with ORFs (tab-separated
      values file with 'ORF', 'Length_log10', 'uATGs', 'FE_atg',
      'FE_cap', 'utr', 'utr_gc', 'polyA' columns) (optional)
    * 't_rna_file': tRNA estimates file (tab-separated values file
      with 'AA', 'Codon', 'tRNA', 'tAI', 'Microarray', 'RNA.seq'
      columns)  (optional)
    * While both 'codon_positions_file' and 't_rna_file' are optional
      either both must be specified or neither must be specified.

    Statistics and figure generation parameters:

    * 'buffer': Length of flanking region around the CDS (default 250)
    * 'count_reads': Scan input, temporary and output files and
      produce counts of reads in each FASTQ, SAM, and BAM file
      processed? (default 'TRUE')
    * 'count_threshold': Remove genes with a read count below this
      threshold, when generating statistics and figures (default 1)
    * 'dataset': Human-readable name of the dataset (default
       'dataset')
    * 'output_metagene_normalized_profile': Calculate position-specific nucleotide
      freqeuency? (default 'TRUE')
    * 'feature': Feature type (default 'CDS')
    * 'is_riboviz_gff': Does the GFF file contain 3 elements per gene
      - UTR5, CDS, and UTR3? (default 'TRUE'). Used by 'bam_to_h5.R'
      only.
    * 'make_bedgraph': Output bedgraph data files in addition to H5
      files? (default 'TRUE')
    * 'max_read_length': Maximum read length in H5 output (default 50)
    * 'min_read_length': Minimum read length in H5 output (default 10)
    * 'primary_id': Primary gene IDs to access the data (YAL001C,
      YAL003W, etc.) (default 'Name')
    * 'rpf': Is the dataset an RPF or mRNA dataset? (default 'TRUE')
    * 'secondary_id': Secondary gene IDs to access the data (COX1,
      EFB1, etc. or 'NULL') (default 'NULL')
    * 'stop_in_feature': Are stop codons part of the feature
      annotations in GFF? Used by 'bam_to_h5.R' only (and only
      if 'is_riboviz_gff' is 'FALSE') (default 'FALSE').

    Visualization parameters:

    * 'run_static_html': run static html visualization per sample? (default 'TRUE')
    * 'output_pdfs': generate .pdfs for sample-related plots (default 'TRUE')

    General:

    * 'validate_only': Validate configuration, check that mandatory
      parameters have been provided and that input files exist, then
      exit without running the workflow? (default 'FALSE')
    * 'publish_index_tmp': Publish/copy index and temporary files to
      'dir_index' and 'dir_tmp'. If 'FALSE' then only symbolic links
      to these files in the Nextflow 'work/' directory are
      created in 'dir_index' and 'dir_tmp' (default 'FALSE')
    * 'skip_inputs': When validating configuration (see
      'validate_only' above) skip checks for existence of ribosome
      profiling data files ('fq_files', 'multiplexed_fq_files',
      'sample_sheet')? (default 'FALSE')
    * 'num_processes': Number of processes to parallelize over, used
      by specific steps in the workflow (default 1)
    * 'samsort_memory': Memory to give to 'samtools sort' (
      default '768M', 'samtools sort' built-in default,
      see http://www.htslib.org/doc/samtools-sort.html)

    Environment variable and configuration tokens:

    * The following configuration parameters take values that are
      absolute or relative paths to files or directories:
      - 'asite_disp_length_file'
      - 'codon_positions_file'
      - 'dir_in'
      - 'dir_index'
      - 'dir_tmp'
      - 'dir_out'
      - 'features_file'
      - 'orf_fasta_file'
      - 'orf_gff_file'
      - 'rrna_fasta_file'
      - 't_rna_file'
    * To give you flexibility in how and where you locate these input
      files and directories, the values for these paths can each
      include, as a prefix, one of the following three tokens:
      - '\${RIBOVIZ_SAMPLES}'
      - '\${RIBOVIZ_ORGANISMS}'
      - '\${RIBOVIZ_DATA}'
    * At runtime, these tokens will be replaced with the names of the
      corresponding environment variables. The environment variables
      are:
      - 'RIBOVIZ_SAMPLES'
      - 'RIBOVIZ_ORGANISMS'
      - 'RIBOVIZ_DATA'
    * If a token is present in a value but the corresponding
      environment variable is undefined, then the path '.' is
      substituted.
    * Which if any token you use in each of the configuration
      parameters is entirely up to you. No checks are made to see
      which specific token is used with which configuration
      parameter.
    * The following configuration parameters also take values that are
      relative paths to files or directories, but the use of tokens in
      their values is *not supported*, as these paths are assumed to
      be relative to paths defined by the configuration parameters
      stated above:
      -'fq_files',' relative to 'dir_in'.
      - 'multiplex_fq_files', relative to 'dir_in'.
      - 'orf_index_prefix', relative to 'dir_index'.
      - 'sample_sheet', relative to 'dir_in'.
      - 'rrna_index_prefix', relative to 'dir_index'.
    """.stripIndent()
    print(help)
}

/**
 * Customise string. Given a string and mapping from tokens to
 * substrings, iterate through the tokens and when a match is
 * found replace the token with substring in the original string,
 * returning the new string. Only the first matching token
 * is replaced.
 * @param string String
 * @param token_replacements Mapping from tokens to substrings
 */
def replace_tokens(string, token_replacements)
{
    for (token_replace in token_replacements)
    {
        if (string.indexOf(token_replace.key) >= 0) {
            return string.replace(token_replace.key, token_replace.value)
        }
    }
    return string
}

// Help message implementation, following
// https://github.com/nf-core/rnaseq/blob/master/main.nf (MIT License)
params.help = false
if (params.help) {
    helpMessage()
    exit 0
}

/*
Initialise and validate configuration.

First initialise optional variables to default values.
This avoids "WARN: Access to undefined parameter '<PARAM>'" errors.
Note that parameter values in the config file override these defaults.
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
params.output_metagene_normalized_profile = true
params.hisat2_orf_params = "-k 2 --no-spliced-alignment --rna-strandness F --no-unal"
params.extract_umis = false
params.trim_5p_mismatches = true
params.feature = "CDS"
params.fq_files = [:]
params.group_umis = false
params.dedup_stats = true
params.is_riboviz_gff = true
params.make_bedgraph = true
params.max_read_length = 50
params.min_read_length = 10
params.multiplex_fq_files = []
params.num_processes = 1
params.output_pdfs = true
params.publish_index_tmp = false
params.primary_id = "Name"
params.rpf = true
params.run_static_html = true
params.secondary_id = null
params.stop_in_feature = false
params.samsort_memory = null
params.validate_only = false
params.skip_inputs = false

if (params.publish_index_tmp)
{
     params.publish_index_tmp_type = 'copy'
}
else
{
    params.publish_index_tmp_type = 'symlink'
}

if (params.validate_only) {
    println("Validating configuration only")
}
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
    exit 1, "CDS flanking region length (buffer) must be >= 0"
}
if (params.count_threshold < 0) {
    exit 1, "Read count threshold (count_threshold) must be >= 0"
}
if (params.num_processes < 1) {
    exit 1, "Number of processes (num_processes) must be >= 1"
}
if (params.min_read_length < 1) {
    exit 1, "Minimum read length in H5 output (min_read_length) must be >= 1"
}
if (params.max_read_length < 1) {
    exit 1, "Maximum read length in H5 output (max_read_length) must be >= 1"
}
if (params.max_read_length < params.min_read_length) {
    exit 1, "Maximum read length in H5 output (max_read_length) must be >= minimum read length (min_read_length)"
}
//if (! params.secondary_id) {
//    secondary_id = null
//} else {
//    secondary_id = params.secondary_id
//}
if (params.dedup_umis) {
    if (! params.extract_umis) {
        println("Warning: deduplication was requested (dedup_umi: TRUE) but UMI extraction was not (extract_umis: FALSE)")
    }
}
if (params.extract_umis) {
    if (! params.containsKey('umi_regexp')) {
        exit 1, "Undefined barcode/UMI regular expression (umi_regexp) when UMI extraction is requested (extract_umis: TRUE)"
    }
}

/*
Get and validate environment variables
*/

samples_dir = System.getenv('RIBOVIZ_SAMPLES') ?: "."
println("samples_dir: " + samples_dir)
if (! file(samples_dir).exists())
{
    exit 1, "No such RIBOVIZ_SAMPLES directory: $samples_dir"
}
organisms_dir = System.getenv('RIBOVIZ_ORGANISMS') ?: "."
println("organisms_dir: " + organisms_dir)
if (! file(organisms_dir).exists())
{
    exit 1, "No such RIBOVIZ_ORGANISMS directory: $organisms_dir"
}
data_dir = System.getenv('RIBOVIZ_DATA') ?: "."
println("data_dir: " + data_dir)
if (! file(data_dir).exists())
{
    exit 1, "No such RIBOVIZ_DATA directory: $data_dir"
}
// Create mapping from environment variable tokens to directories.
riboviz_env_paths = [:]
riboviz_env_paths['\${RIBOVIZ_SAMPLES}'] = samples_dir
riboviz_env_paths['\${RIBOVIZ_ORGANISMS}'] = organisms_dir
riboviz_env_paths['\${RIBOVIZ_DATA}'] = data_dir

/*
Apply environment variables to paths
*/

if (! params.containsKey('dir_in')) {
    exit 1, "Undefined input directory (dir_in)"
}
params.dir_in_env = replace_tokens(params.dir_in, riboviz_env_paths)
params.dir_index_env = replace_tokens(params.dir_index, riboviz_env_paths)
params.dir_out_env = replace_tokens(params.dir_out, riboviz_env_paths)
params.dir_tmp_env = replace_tokens(params.dir_tmp, riboviz_env_paths)

/*
Validate input files.
*/

sample_id_fq = [:]
multiplex_id_fq = [:]
multiplex_sample_sheet_tsv = Channel.empty()
is_multiplexed = false
if (params.validate_only && params.skip_inputs) {
    println("Skipping checks for existence of of ribosome profiling input files (params.dir_in_env|fq_files|multiplex_fq_files|sample_sheet)")
}
if ((! params.validate_only) || (! params.skip_inputs)) {
    if (! file(params.dir_in_env).exists())
    {
        exit 1, "No such directory (dir_in): $params.dir_in_env"
    }
}
if ((! params.fq_files) && (! params.multiplex_fq_files)) {
    exit 1, "No sample files (fq_files) or multiplexed files (multiplex_fq_files) are defined"
} else if (params.fq_files && params.multiplex_fq_files) {
    exit 1, "Both sample files (fq_files) and multiplexed files (multiplex_fq_files) are defined - only one or the other should be defined"
} else if (params.fq_files) {
    if ((! params.validate_only) || (! params.skip_inputs)) {
        // Filter 'params.fq_files' down to those samples that exist.
        for (entry in params.fq_files) {
            sample_fq = file("${params.dir_in_env}/${entry.value}")
            if (sample_fq.exists()) {
                sample_id_fq[entry.key] = sample_fq
            } else {
                println("No such sample file ($entry.key): $entry.value")
            }
        }
        if (! sample_id_fq) {
            exit 1, "None of the defined sample files (fq_files) exist"
        }
    }
} else {
    if ((! params.validate_only) || (! params.skip_inputs)) {
        // Filter 'params.multiplex_fq_files' down to those files that exist.
        for (entry in params.multiplex_fq_files) {
            multiplex_fq = file("${params.dir_in_env}/${entry}")
            if (multiplex_fq.exists()) {
                // Use file base name as key, ensuring that if file
                // has extension '.fastq.gz' or '.fq.gz' then both
                // extensions are removed from the name.
                multiplex_id = multiplex_fq.baseName
                if (multiplex_id.endsWith(".fastq")) {
                    multiplex_id = multiplex_id - '.fastq'
                } else if (multiplex_id.endsWith(".fq")) {
                    multiplex_id = multiplex_id - '.fq'
                }
                multiplex_id_fq[multiplex_id] = multiplex_fq
            } else {
                println("No such multiplexed file: $entry")
            }
        }
        if (! multiplex_id_fq) {
            exit 1, "None of the defined multiplexed files (multiplex_fq_files) exist"
        }
    }
    if (! params.containsKey('sample_sheet')) {
        exit 1, "Undefined sample sheet (sample_sheet)"
    }
    sample_sheet = file("${params.dir_in_env}/${params.sample_sheet}")
    if ((! params.validate_only) || (! params.skip_inputs)) {
        if (! sample_sheet.exists()) {
            exit 1, "No such sample sheet (sample_sheet): ${sample_sheet}"
        }
        multiplex_sample_sheet_tsv = Channel.fromPath(sample_sheet,
                                                      checkIfExists: true)
    }
    is_multiplexed = true
}

// Create YAML fragment including 'params.fq_files' and
// 'params.multiplex_fq_files' to serve as a configuration file for
// riboviz.tools.count_reads. There is no way to get the location
// of the YAML configuration file itself (i.e. the value of
// '-param-file') from within Nextflow so this is a workaround.
Map ribosome_fqs = [:]
ribosome_fqs.fq_files = params.fq_files
ribosome_fqs.multiplex_fq_files = params.multiplex_fq_files
ribosome_fqs_yaml = new Yaml().dump(ribosome_fqs)

// Non-sample-specific input files.
if (! params.build_indices) {
    rrna_index_prefix = file("${params.dir_index_env}/${params.rrna_index_prefix}.*.ht2")
    if (! rrna_index_prefix) {
        exit 1, "No such rRNA index files (rrna_index_prefix): ${params.dir_index_env}/${params.rrna_index_prefix}.*.ht2"
    }
    pre_built_rrna_index_ht2 = Channel
        .fromPath(rrna_index_prefix, checkIfExists: true)
        .collect()
    orf_index_prefix = file("${params.dir_index_env}/${params.orf_index_prefix}.*.ht2")
    if (! orf_index_prefix) {
        exit 1, "No such ORF index files (orf_index_prefix): ${params.dir_index_env}/${params.orf_index_prefix}.*.ht2"
    }
    pre_built_orf_index_ht2 = Channel
        .fromPath("${params.dir_index_env}/${params.orf_index_prefix}.*.ht2",
                  checkIfExists: true)
        .collect()
} else {
    pre_built_rrna_index_ht2 = Channel.empty()
    pre_built_orf_index_ht2 = Channel.empty()
}

if (! params.containsKey('rrna_fasta_file')) {
    exit 1, "Undefined rRNA FASTA file (rrna_fasta_file)"
}
rrna_fasta_file = file(replace_tokens(params.rrna_fasta_file,
                                      riboviz_env_paths))
if (! rrna_fasta_file.exists()) {
    exit 1, "No such file rRNA FASTA file (rrna_fasta_file): ${rrna_fasta_file}"
}
rrna_fasta = Channel.fromPath(rrna_fasta_file, checkIfExists: true)
if (! params.containsKey('orf_fasta_file')) {
    exit 1, "Undefined ORF FASTA file (orf_fasta_file)"
}
orf_fasta_file = file(replace_tokens(params.orf_fasta_file,
                                     riboviz_env_paths))
if (! orf_fasta_file.exists()) {
    exit 1, "No such ORF FASTA file (orf_fasta_file): ${orf_fasta_file}"
}
orf_fasta = Channel.fromPath(orf_fasta_file, checkIfExists: true)
if (! params.containsKey('orf_gff_file')) {
    exit 1, "Undefined ORF GFF file (orf_gff_file)"
}
orf_gff_file = file(replace_tokens(params.orf_gff_file, riboviz_env_paths))
if (! orf_gff_file.exists()) {
    exit 1, "No such ORF GFF file (orf_gff_file): ${orf_gff_file}"
}
orf_gff = Channel.fromPath(orf_gff_file, checkIfExists: true)

// Optional inputs for generate_stats_figs.R.
// If an optional file is not provided then a 'Missing_<PARAM>' file
// (for example 'Missing_features_file') is created within the 'work/'
// directories for the generateStatsFigs process. This symbolically
// links to a non-existent 'Missing_<PARAM>' file in the users current
// directory. This is not an issue since the files will not be passed
// onto generate_stats_figs.R and no attempt is made to use them. They
// are a side-effect of using the Nextflow pattern for optional
// inputs.
// Optional inputs implementation follows pattern
// https://github.com/nextflow-io/patterns/blob/master/optional-input.nf.
if (params.containsKey('t_rna_file') && params.t_rna_file) {
    t_rna_file = file(replace_tokens(params.t_rna_file,
                                     riboviz_env_paths))
    if (! t_rna_file.exists()) {
        exit 1, "No such tRNA estimates file (t_rna_file): ${t_rna_file}"
    }
    t_rna_tsv = Channel.fromPath(t_rna_file, checkIfExists: true)
    params.is_t_rna_file = true
} else {
    t_rna_tsv = file("Missing_t_rna_file")
    params.is_t_rna_file = false
}
if (params.containsKey('codon_positions_file')
    && params.codon_positions_file) {
    codon_positions_file = file(replace_tokens(params.codon_positions_file,
                                               riboviz_env_paths))
    if (! codon_positions_file.exists()) {
        exit 1, "No such codon positions file (codon_positions_file): ${codon_positions_file}"
    }
    codon_positions_rdata = Channel.fromPath(codon_positions_file,
                                             checkIfExists: true)
    params.is_codon_positions_file = true
} else {
    codon_positions_rdata = file("Missing_codon_positions_file")
    params.is_codon_positions_file = false
}
if (params.is_t_rna_file && params.is_codon_positions_file) {
    params.is_t_rna_and_codon_positions_file = true
} else if ((! params.is_t_rna_file) && (! params.is_codon_positions_file)) {
    params.is_t_rna_and_codon_positions_file = false
} else {
    exit 1, "Either both tRNA estimates (t_rna_file) and codon positions (codon_positions_file) must be defined or neither must be defined"
}
if (params.containsKey('features_file') && params.features_file) {
    features_file = file(replace_tokens(params.features_file,
                                        riboviz_env_paths))
    if (! features_file.exists()) {
        exit 1, "No such features file (features_file): ${features_file}"
    }
    features_tsv = Channel.fromPath(features_file, checkIfExists: true)
    params.is_features_file = true
} else {
    features_tsv = file("Missing_features_file")
    params.is_features_file = false
}
if (params.containsKey('asite_disp_length_file')
    && params.asite_disp_length_file) {
    asite_disp_length_file = file(replace_tokens(params.asite_disp_length_file,
                                                 riboviz_env_paths))
    if (! asite_disp_length_file.exists()) {
        exit 1, "No such A-site displacement file (asite_disp_length_file): ${asite_disp_length_file}"
    }
    asite_disp_length_txt = Channel.fromPath(asite_disp_length_file,
                                             checkIfExists: true)
    params.is_asite_disp_length_file = true
} else {
    asite_disp_length_txt = file("Missing_aside_disp_length_file")
    params.is_asite_disp_length_file = false
}

if (params.validate_only) {
    exit 0, "Validated configuration"
}

include { buildIndicesrRNA; buildIndicesORF } from './modules/build_indices'
include { cutAdapters; extractUmis; cutAdaptersMultiplex; extractUmisMultiplex; demultiplex } from './modules/preprocess'
include { hisat2rRNA; hisat2ORF } from './modules/alignment'
include { trim5pMismatches; samViewSort; groupUmisPreDedup; dedupUmis; groupUmisPostDedup; outputBams; makeBedgraphs; bamToH5 } from './modules/postprocess'
include { generateStatsFigs; renameTpms; collateTpms; createVizParamsConfigFile; staticHTML; createInteractiveVizParamsConfigFile; countReads } from './modules/visualization_and_stats'


workflow buildIndices {

    take:
    rrna_fasta
    orf_fasta

    main:
    built_rrna_index_ht2 = buildIndicesrRNA(rrna_fasta)
    built_orf_index_ht2 = buildIndicesORF(orf_fasta)

    emit:
    built_rrna_index_ht2
    built_orf_index_ht2

}

/*
Sample file (fq_files)-specific processes.
*/

workflow preprocessReads{

    take:
    sample_id_fq

    main:
    
    cut_fq = cutAdapters(Channel.from(sample_id_fq.collect{ id, file -> tuple(id, file) }))
    // Route 'cut_fq' channel outputs depending on whether UMIs are to be
    // extracted or not.
    cut_fq.branch {
        umi_fq: params.extract_umis
        non_umi_fq: ! params.extract_umis
    }
    .set { cut_fq_branch }

    umi_extract_fq = extractUmis(cut_fq_branch.umi_fq)
    
    trimmed_fq = cut_fq_branch.non_umi_fq
          .mix(umi_extract_fq)

    emit:
    trimmed_fq
}

workflow preprocessMultiplexedReads{

      take:
      multiplex_id_fq

      main:

      cut_multiplex_fq = cutAdaptersMultiplex(Channel.from(multiplex_id_fq.collect{ id, file -> tuple(id, file) }))

      // Route 'cut_multiplex_fq' channel outputs depending on whether UMIs
      // are to be extracted or not.
      cut_multiplex_fq.branch {
          umi_fq: params.extract_umis
          non_umi_fq: ! params.extract_umis
      }
      .set { cut_multiplex_fq_branch }

      umi_extract_multiplex_fq = extractUmis(cut_multiplex_fq_branch.umi_fq)
      
      // Combine channels for downstream processing. By definition of
      // 'cut_multiplex_fq.branch' only one of the input channels will have
      // content.
      trimmed_multiplex_fq = cut_multiplex_fq_branch.non_umi_fq
          .mix(umi_extract_multiplex_fq)

      // Use '.toString' to prevent changing hashes of
      // 'workflow.projectDir' triggering reexecution of this
      // process if 'nextflow run' is run with '-resume'.
      demultiplex(workflow.projectDir.toString(),trimmed_multiplex_fq, multiplex_sample_sheet_tsv)
      demultiplex_fq = demultiplex.out.demultiplex_fq
      demultiplex_num_reads_tsv = demultiplex.out.demultiplex_num_reads_tsv
      // 'demultiplex_fq' outputs a single list with all the output
      // files. Extract sample IDs from file basenames, filter out
      // 'Unassigned' and output tuples of sample IDs and file names as
      // separate items onto a new channel.
      demultiplex_samples_fq = demultiplex_fq
            .flatten()
            // Use file basename as sample ID.
            .map { [it.baseName, it] }
            // If file was '.fastq|fq.gz' then basename will include
            // '.fq|fastq' so strip that off too.
            .map { n, f -> [n.endsWith(".fq") ? n - ".fq" : n, f] }
            .map { n, f -> [n.endsWith(".fastq") ? n - ".fastq" : n, f] }
            .filter { n, f -> n != "Unassigned" }
            //.into { report_demultiplex_samples_fq; demultiplex_samples_fq }

      demultiplex_sample_ids = demultiplex_samples_fq
          .map { n, f -> n }
          .toList() // [] if none
          .view { "Demultiplexed samples: ${it}"}
          // Wrap list in list, so 'merge' below doesn't append lists
          .map { it -> [it] }

      multiplex_sample_sheet_ids = multiplex_sample_sheet_tsv
          // Extract original sample IDs from sample sheet
          .splitCsv(header: true, sep: '\t')
          .map { row -> row.SampleID } // No output if no 'SampleID' column
          .toList() // [] if no 'SampleID' column
          // Wrap list in list, so 'merge' below doesn't append lists
          .map { it -> [it] }

      multiplex_sample_sheet_ids
          .merge(demultiplex_sample_ids)
          .map { a, b -> a - b}
          .view { "Non-demultiplexed samples: ${it}" }


      /*
      Sample-specific processes.

      Common to both sample files (fq_files) and demultiplexed files.
      */

      // Combine channels for downstream processing. By definition of
      // upstream conditions and processes, only one of the channels
      // will have content.
      
      //umi_extract_multiplex_fq | view
      
      //trimmed_fq = cut_multiplex_fq_branch.non_umi_fq
      //    .mix(umi_extract_multiplex_fq)
      //    .mix(demultiplex_samples_fq)
          
      trimmed_fq = demultiplex_samples_fq
      

    emit:
    trimmed_fq

}


workflow postProcessMappedReads{
  
    take:
    trim_5p_mismatches
    orf_gff
  
    main:
    // Route 'trim_5p_branch' channel outputs depending on whether mismatched
      // 5' base are to be trimmed or not
      trim_5p_mismatches.branch {
          trim_5p_fq: params.trim_5p_mismatches
          non_trim_5p_fq: ! params.trim_5p_mismatches
      }
      .set { trim_5p_branch }
  
  
      // Use '.toString' to prevent changing hashes of
      // 'workflow.projectDir' triggering reexecution of this
      // process if 'nextflow run' is run with '-resume'.
      trim5pMismatches(workflow.projectDir.toString(),trim_5p_branch.trim_5p_fq)
      trim_orf_map_sam = trim5pMismatches.out.trim_orf_map_sam


      // Combine channels for downstream processing. By definition of
      // upstream conditions and processes, only one of the channels
      // will have content.
      trimmed_5p_fq = trim_5p_branch.non_trim_5p_fq
          .mix(trim_orf_map_sam)
  
      orf_map_bam = samViewSort(trimmed_5p_fq)
      // Route "orf_map_bam" channel outputs depending on whether UMIs are
      // to be deduplicated or not.
      orf_map_bam.branch {
          dedup_bam: params.dedup_umis
          non_dedup_bam: ! params.dedup_umis
      }
      .set { orf_map_bam_branch }
      pre_output_bam = Channel.empty()
      if (params.dedup_umis && params.group_umis)
      {
        groupUmisPreDedup(orf_map_bam_branch.dedup_bam)
        dedup_bam = dedupUmis(orf_map_bam_branch.dedup_bam)
        groupUmisPostDedup(dedup_bam.dedup_bam)
        pre_output_bam = dedup_bam.dedup_bam
      } else {
        pre_output_bam = orf_map_bam_branch.non_dedup_bam
      }
      output_bam = outputBams(pre_output_bam)
      if (params.make_bedgraph)
      {
        makeBedgraphs(output_bam)
      }
  
      h5s = bamToH5(output_bam,orf_gff)

    emit:
    h5s

}


workflow visualizeResults{
  take:
  h5s
  orf_fasta
  orf_gff
  t_rna_tsv
  codon_positions_rdata
  features_tsv
  asite_disp_length_txt


  main:
  
  missing_options = h5s.map{ sample_id, sample_h5, sample_h5_star -> sample_id}.combine(Channel.fromPath("."))
  
  generateStatsFigs(h5s,orf_fasta,orf_gff,t_rna_tsv,codon_positions_rdata,features_tsv,asite_disp_length_txt)
  
  if (params.is_asite_disp_length_file) {
        read_frame_per_orf_filtered_tsv = generateStatsFigs.out.read_frame_per_orf_filtered_tsv
  } else {
        read_frame_per_orf_filtered_tsv = missing_options
  }
  if (params.is_t_rna_and_codon_positions_file) {
        normalized_density_apesites_per_codon_long_tsv = generateStatsFigs.out.normalized_density_apesites_per_codon_long_tsv
  } else {
        normalized_density_apesites_per_codon_long_tsv = missing_options
  }
  if (params.is_features_file) {
        ORF_TPMs_vs_features_tsv = generateStatsFigs.out.ORF_TPMs_vs_features_tsv
  } else {
        ORF_TPMs_vs_features_tsv = missing_options
  }
  
  
  
  // Join outputs from generateStatsFigs for staticHTML.
  // Join is done on first value of each tuple i.e. sample ID.

  generate_stats_figs_static_html =
      generateStatsFigs.out.metagene_start_stop_read_counts_tsv
      .join(generateStatsFigs.out.metagene_position_length_counts_5start_tsv, remainder: true)
      .join(generateStatsFigs.out.read_counts_by_length_tsv, remainder: true)
      .join(generateStatsFigs.out.metagene_normalized_profile_start_stop_tsv, remainder: true)
      .join(read_frame_per_orf_filtered_tsv, remainder: true)
      .join(ORF_TPMs_vs_features_tsv, remainder: true)
      .join(normalized_density_apesites_per_codon_long_tsv, remainder: true)
  
  
  generateStatsFigs.out.finished_sample_id
      .ifEmpty { exit 1, "No sample was processed successfully" }
      .view { "Finished processing sample: ${it}" }

  renameTpms(generateStatsFigs.out.orf_tpms_and_counts_tsv)
  collateTpms(renameTpms.out.tpms_sample_id.collect(),renameTpms.out.tpms_sample_tsv.collect())

  count_reads_sample_ids = collateTpms.out.collate_tpms_sample_ids.collect()

  emit:
  generate_stats_figs_static_html
  count_reads_sample_ids
}
  

workflow generateHTML {
  take:
  generate_stats_figs_static_html
  viz_params_yaml
  interactive_viz_params_yaml

  main:

  viz_params_config_file_yaml = createVizParamsConfigFile(viz_params_yaml)
  staticHTML(viz_params_config_file_yaml,generate_stats_figs_static_html)

  createInteractiveVizParamsConfigFile(interactive_viz_params_yaml)
  static_html_sample_ids = staticHTML.out.static_html_sample_ids.collect()
  
  // Create handler for finished_viz_sample_id channel, output by
  // staticHTML, only if run_static_html is true i.e. if staticHTML
  // executes.
  if (params.run_static_html) {
    staticHTML.out.finished_viz_sample_id
        .ifEmpty { exit 1, "No sample was visualised successfully" }
         .view { "Finished visualising sample: ${it}" }
  }
  
  emit:
  static_html_sample_ids

}


workflow finalReadCount {
  take:
  ribosome_fqs_yaml
  count_reads_sample_ids

  main:
  // Use '.toString' to prevent changing hashes of
  // 'workflow.projectDir' triggering reexecution of this
  // process if 'nextflow run' is run with '-resume'.
  countReads(workflow.projectDir.toString(),ribosome_fqs_yaml,count_reads_sample_ids)

}

workflow {
   
    if (params.build_indices)
    {
      buildIndices(rrna_fasta,orf_fasta)
      rrna_index_ht2 = buildIndices.out.built_rrna_index_ht2
      orf_index_ht2 = buildIndices.out.built_orf_index_ht2
    } 
    else
    {
      rrna_index_ht2 = pre_built_rrna_index_ht2
      orf_index_ht2 = pre_built_orf_index_ht2
    }
    if (!is_multiplexed)
    {
      trimmed_fq = preprocessReads(sample_id_fq)
    } else {
      trimmed_fq = preprocessMultiplexedReads(multiplex_id_fq)
    }

    hisat2rRNA(trimmed_fq,rrna_index_ht2)
    hisat2ORF(hisat2rRNA.out.non_rrna_fq,orf_index_ht2)
    trim_5p_mismatches = hisat2ORF.out.trim_5p_mismatches

    h5s = postProcessMappedReads(trim_5p_mismatches,orf_gff)
    visualizeResults(h5s,orf_fasta,orf_gff,t_rna_tsv,codon_positions_rdata,features_tsv,asite_disp_length_txt)

    Map viz_params = [:]
    if (params.is_asite_disp_length_file) {
        viz_params.asite_disp_length_file = asite_disp_length_file.toString()
    }
    if (params.is_codon_positions_file) {
        viz_params.codon_positions_file = codon_positions_file.toString()
    }
    if (params.is_features_file) {
        viz_params.features_file = features_file.toString()
    }
    if (params.is_t_rna_file) {
        viz_params.t_rna_file = t_rna_file.toString()
    }

   viz_params_yaml = new Yaml().dump(viz_params)

    // collect only parameters needed for interactive visualization (riboviz/#275)
    // NOTE: fq_files, dataset & sample_sheet don't use environment tokens, are relative to dir_in
    // however dir_in, dir_out and features_file MAY use environment tokens but these are handled above in the script
    Map interactive_viz_params = [:]
    interactive_viz_params.dir_in = params.dir_in_env
    interactive_viz_params.dir_out = params.dir_out_env
    interactive_viz_params.dataset = params.dataset
    interactive_viz_params.fq_files = params.fq_files
    interactive_viz_params.sample_sheet = params.sample_sheet
    if (params.is_features_file) {
        interactive_viz_params.features_file = features_file.toString()
    }
    interactive_viz_params_yaml = new Yaml().dump(interactive_viz_params)

    if (params.run_static_html) {
      generateHTML(visualizeResults.out.generate_stats_figs_static_html,viz_params_yaml,interactive_viz_params_yaml)
      count_reads_sample_ids = generateHTML.out.static_html_sample_ids
    } else {
      count_reads_sample_ids = visualizeResults.out.count_reads_sample_ids
    }
    
    if (params.count_reads)
    {
      finalReadCount(ribosome_fqs_yaml,count_reads_sample_ids)
    }

    
}



workflow.onComplete {
    println "Workflow finished! (${workflow.success ? 'OK' : 'failed'})"
}
