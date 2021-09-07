"""
Configuration parameter names.
"""

INPUT_DIR = "dir_in"
""" Input directory. """
INDEX_DIR = "dir_index"
""" Index files directory. """
TMP_DIR = "dir_tmp"
""" Intermediate files directory. """
OUTPUT_DIR = "dir_out"
""" Output files directory. """

OUTPUT_PDFS = "output_pdfs"
""" Generate .pdfs for sample-related plots. """

ORF_FASTA_FILE = "orf_fasta_file"
""" ORF file to align to. """
ORF_GFF_FILE = "orf_gff_file"
""" GFF2/GFF3 file for ORFs. """
RRNA_FASTA_FILE = "rrna_fasta_file"
""" rRNA file to avoid aligning to. """
CODON_POSITIONS_FILE = "codon_positions_file"
""" Codon positions in each gene. """
FEATURES_FILE = "features_file"
""" Features to correlate with ORFs. """
T_RNA_FILE = "t_rna_file"
""" tRNA estimates. """
ASITE_DISP_LENGTH_FILE = "asite_disp_length_file"
""" Fixed A-site positions by read length. """

FQ_FILES = "fq_files"
" fastq files to be processed. """
MULTIPLEX_FQ_FILES = "multiplex_fq_files"
" Multiplexed fastq files to be processed. """

TRIM_5P_MISMATCHES = "trim_5p_mismatches"
""" Trim mismatched 5' base? """
RUN_STATIC_HTML = "run_static_html"
""" Create static html visualization per sample? """

SAMPLE_SHEET = "sample_sheet"
""" Sample sheet. """
DEDUP_UMIS = "dedup_umis"
""" Extract UMIs and deduplicate reads? """
GROUP_UMIS = "group_umis"
""" Summarise UMI groups before and after deduplication? """
DEDUP_STATS = "dedup_stats"
""" Output UMI deduplication statistics? """
EXTRACT_UMIS = "extract_umis"
""" Extract UMIs? """
UMI_REGEXP = "umi_regexp"
"""
UMI-tools-compliant regular expression to extract UMIs and
barcodes.
"""

BUILD_INDICES = "build_indices"
""" Build indices for aligner flag. """
ORF_INDEX_PREFIX = "orf_index_prefix"
""" rRNA index file name prefix. """
RRNA_INDEX_PREFIX = "rrna_index_prefix"
""" ORF index file name prefix. """

ADAPTERS = "adapters"
""" Illumina sequencing adapter to remove. """
MAKE_BEDGRAPH = "make_bedgraph"
""" Output bedgraph files. """
BUFFER = "buffer"
""" Length of flanking region around the CDS. """
COUNT_THRESHOLD = "count_threshold"
"""
Remove genes with a read count below this threshold, when generating
statistics and figures.
"""
DATASET = "dataset"
""" Dataset name. """
OUTPUT_METAGENE_NORMALIZED_PROFILE = "output_metagene_normalized_profile"
""" Calculate position-specific nucleotide frequency? """
FEATURE = "feature"
""" Feature type """
MIN_READ_LENGTH = "min_read_length"
""" Minimum read length in H5 output. """
MAX_READ_LENGTH = "max_read_length"
""" Maximum read length in H5 output. """
PRIMARY_ID = "primary_id"
""" Primary gene IDs. """
SECONDARY_ID = "secondary_id"
""" Secondary gene IDs. """
IS_RIBOVIZ_GFF = "is_riboviz_gff"
"""
Does the GFF file contain 3 elements per gene - UTR5, CDS, and UTR3?
"""
RPF = "rpf"
""" Is the dataset an RPF or mRNA dataset? """
STOP_IN_FEATURE = "stop_in_feature"
""" Are stop codons part of the feature annotations in GFF? """
COUNT_READS = "count_reads"
"""
Scan input, temporary and output files and produce counts of reads in
each FASTQ, SAM and BAM file processed?
"""

NUM_PROCESSES = "num_processes"
""" Number of processes to parallelize over. """
SAMSORT_MEMORY = "samsort_memory"
""" Memory to give to 'samtools sort'. """

VALIDATE_ONLY = "validate_only"
""" Validate configuration only? """
PUBLISH_INDEX_TMP = "publish_index_tmp"
""" Publish index and temporary files? """
SKIP_INPUTS = "skip_inputs"
"""
When validating configuration skip checks for existence of ribosome
profiling data files.
"""

DEFAULT_FOLDER_VALUES = {
    INDEX_DIR: "index",
    TMP_DIR: "tmp",
    OUTPUT_DIR: "output"
}
"""
Default values for folder parameters.
Consistent with ``prep_riboviz.nf``.
"""

DEFAULT_CONDITIONS = {
    BUILD_INDICES: True,
    TRIM_5P_MISMATCHES: True,
    EXTRACT_UMIS: False,
    DEDUP_UMIS: False,
    DEDUP_STATS: True,
    GROUP_UMIS: False,
    MAKE_BEDGRAPH: True,
    RUN_STATIC_HTML: True,
    COUNT_READS: True,
    OUTPUT_PDFS: True,
    OUTPUT_METAGENE_NORMALIZED_PROFILE: True,
    PUBLISH_INDEX_TMP: False
}
"""
Default values for conditional boolean parameters that affect which
steps of the workflow are invoked and the files output.
Consistent with  ``prep_riboviz.nf``.
"""

R_LIBS = "r_libs"
""" R libraries directory (job submission) (command-line only). """
CONFIG_FILE = "config_file"
"""
RiboViz YAML configuration file (job submission) (command-line only).
"""
NEXTFLOW_DAG_FILE = "nextflow_dag_file"
"""
Nextflow DAG file (job submission).
"""
NEXTFLOW_REPORT_FILE = "nextflow_report_file"
"""
Nextflow report file (job submission).
"""
NEXTFLOW_TIMELINE_FILE = "nextflow_timeline_file"
"""
Nextflow timeline file (job submission).
"""
NEXTFLOW_TRACE_FILE = "nextflow_trace_file"
"""
Nextflow trace file (job submission).
"""
NEXTFLOW_WORK_DIR = "nextflow_work_dir"
"""
Nextflow work directory (job submission).
"""
NEXTFLOW_RESUME = "nextflow_resume"
"""
Resume Nextflow workflow (job submission).
"""
JOB_NAME = "job_name"
""" Name of batch job (job submission). """
JOB_RUNTIME = "job_runtime"
"""
Maximum runtime for batch job (job submission).
"""
JOB_MEMORY = "job_memory"
""" Requested memory for batch job (job submission). """
JOB_NUM_CPUS = "job_num_cpus"
"""
Requested number of CPUs for batch job (job submission).
"""
JOB_PARALLEL_ENV = "job_parallel_env"
""" Requested parallel environment for batch job (job submission). """
JOB_EMAIL = "job_email"
""" E-mail address for batch job events (job submission). """
JOB_EMAIL_EVENTS = "job_email_events"
"""
Events triggering emails about batch job (job submission).
"""

DEFAULT_JOB_CONFIG = {
    JOB_EMAIL: None,
    JOB_EMAIL_EVENTS: "beas",
    JOB_MEMORY: "8G",
    JOB_NAME: "riboviz",
    JOB_NUM_CPUS: 4,
    JOB_RUNTIME: "48:00:00",
    JOB_PARALLEL_ENV: "mpi",
    NEXTFLOW_DAG_FILE: "nextflow-dag.html",
    NEXTFLOW_REPORT_FILE: "nextflow-report.html",
    NEXTFLOW_RESUME: "false",
    NEXTFLOW_TIMELINE_FILE: "nextflow-timeline.html",
    NEXTFLOW_TRACE_FILE: "nextflow-trace.tsv",
    NEXTFLOW_WORK_DIR: "work",
    VALIDATE_ONLY: False
}
""" Default values for job configuration parameters. """
JOB_CONFIG_TYPE = {
    JOB_EMAIL: str,
    JOB_EMAIL_EVENTS: str,
    JOB_MEMORY: str,
    JOB_NAME: str,
    JOB_PARALLEL_ENV: str,
    JOB_RUNTIME: str,
    JOB_NUM_CPUS: int,
    NEXTFLOW_DAG_FILE: str,
    NEXTFLOW_REPORT_FILE: str,
    NEXTFLOW_RESUME: bool,
    NEXTFLOW_TIMELINE_FILE: str,
    NEXTFLOW_TRACE_FILE: str,
    NEXTFLOW_WORK_DIR: str,
    VALIDATE_ONLY: bool
}
""" Types of job configuration parameters. """

ENV_RIBOVIZ_SAMPLES = "RIBOVIZ_SAMPLES"
""" Samples directory environment variable name. """
ENV_RIBOVIZ_ORGANISMS = "RIBOVIZ_ORGANISMS"
""" Organisms directory environment variable name. """
ENV_RIBOVIZ_DATA = "RIBOVIZ_DATA"
""" Data directory environment variable name. """
ENV_DIRS = [
    ENV_RIBOVIZ_SAMPLES,
    ENV_RIBOVIZ_ORGANISMS,
    ENV_RIBOVIZ_DATA
]
"""
Environment variables which can be used to define directories.
"""
ENV_INPUT_PARAMS = [ASITE_DISP_LENGTH_FILE,
                    CODON_POSITIONS_FILE,
                    FEATURES_FILE,
                    INPUT_DIR,
                    ORF_FASTA_FILE,
                    ORF_GFF_FILE,
                    RRNA_FASTA_FILE,
                    T_RNA_FILE]
"""
Names of input parameters whose values can include RiboViz environment
variables.
"""
ENV_OUTPUT_PARAMS = [INDEX_DIR, OUTPUT_DIR, TMP_DIR]
"""
Names of output parameters whose values can include RiboViz environment
variables.
"""
ENV_PARAMS = ENV_INPUT_PARAMS + ENV_OUTPUT_PARAMS
"""
Names of parameters whose values can include RiboViz environment
variables.
"""
