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

CMD_FILE = "cmd_file"
""" Bash commands file name. """
LOGS_DIR = "dir_logs"
""" Log files directory. """

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
DO_POS_SP_NT_FREQ = "do_pos_sp_nt_freq"
""" Calculate position-specific nucleotide frequency? """
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
STOP_IN_CDS = "stop_in_cds"
""" Are stop codons part of the CDS annotations in GFF? """
COUNT_READS = "count_reads"
"""
Scan input, temporary and output files and produce counts of reads in
each FASTQ, SAM and BAM file processed?
"""

NUM_PROCESSES = "num_processes"
""" Number of processes to parallelize over. """
IS_TEST_RUN = "is_test_run"
""" Is this a test run? (unused). """
ALIGNER = "aligner"
""" Short read aligner to use (unused). """

VALIDATE_ONLY = "validate_only"
""" Validate configuration only? (Nextflow workflow only). """
PUBLISH_INDEX_TMP = "publish_index_tmp"
"""
Publish index and temporary files? (Nextflow workflow only).
"""
SKIP_INPUTS = "skip_inputs"
"""
When validating configuration skip checks for existence of ribosome
profiling data files.
"""

R_LIBS = "r_libs"
""" R libraries directory (job submission) (command-line only). """
CONFIG_FILE = "config_file"
"""
RiboViz YAML configuration file (job submission) (command-line only).
"""
NEXTFLOW_REPORT_FILE = "nextflow_report_file"
""" Nextflow report file (Nextflow workflow only)
(job submission).
"""
NEXTFLOW_WORK_DIR = "nextflow_work_dir"
""" Nextflow work directory (Nextflow workflow only)
(job submission).
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
JOB_EMAIL = "job_email"
""" E-mail address for batch job events (job submission). """
JOB_EMAIL_EVENTS = "job_email_events"
"""
Events triggering emails about batch job (job submission).
"""

DEFAULT_JOB_CONFIG = {
    JOB_NAME: "riboviz",
    JOB_RUNTIME: "48:00:00",
    JOB_MEMORY: "8GB",
    JOB_NUM_CPUS: 4,
    JOB_EMAIL_EVENTS: "beas",
    JOB_EMAIL: None,
    NEXTFLOW_WORK_DIR: "work",
    NEXTFLOW_REPORT_FILE: "nextflow-report.html",
    VALIDATE_ONLY: False
}
""" Default values for job configuration parameters. """
JOB_CONFIG_TYPE = {
    JOB_NAME: str,
    JOB_RUNTIME: str,
    JOB_MEMORY: str,
    JOB_NUM_CPUS: int,
    JOB_EMAIL_EVENTS: str,
    JOB_EMAIL: str,
    NEXTFLOW_WORK_DIR: str,
    NEXTFLOW_REPORT_FILE: str,
    VALIDATE_ONLY: bool
}
""" Types of job configuration parameters. """
