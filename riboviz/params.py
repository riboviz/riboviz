"""
RiboViz configuration parameters.
"""

# Input and output directories.
CMD_FILE = "cmd_file"
INPUT_DIR = "dir_in"
INDEX_DIR = "dir_index"
LOGS_DIR = "dir_logs"
OUTPUT_DIR = "dir_out"
TMP_DIR = "dir_tmp"

# Input files.
FEATURES_FILE = "features_file"
ORF_FASTA_FILE = "orf_fasta_file"
ORF_GFF_FILE = "orf_gff_file"
R_RNA_FASTA_FILE = "rrna_fasta_file"

# Indexing.
BUILD_INDICES = "build_indices"
ORF_INDEX_PREFIX = "orf_index_prefix"
R_RNA_INDEX_PREFIX = "rrna_index_prefix"

# Adapter trimming.
ADAPTERS = "adapters"

# Sample files.
FQ_FILES = "fq_files"

# Demultiplexing and deduplication.
SAMPLE_SHEET = "sample_sheet"
DEDUP_UMIS = "dedup_umis"
GROUP_UMIS = "group_umis"
EXTRACT_UMIS = "extract_umis"
MULTIPLEX_FQ_FILES = "multiplex_fq_files"
UMI_REGEXP = "umi_regexp"

# Bedgraphs.
MAKE_BEDGRAPH = "make_bedgraph"

# Bam to H5 and statistics generation.
BUFFER = "buffer"
CODON_POSITIONS_FILE = "codon_positions_file"
COUNT_THRESHOLD = "count_threshold"
DATASET = "dataset"
DO_POS_SP_NT_FREQ = "do_pos_sp_nt_freq"
MAX_READ_LENGTH = "max_read_length"
MIN_READ_LENGTH = "min_read_length"
PRIMARY_ID = "primary_id"
IS_RIBOVIZ_GFF = "is_riboviz_gff"
RPF = "rpf"
SECONDARY_ID = "secondary_id"
STOP_IN_CDS = "stop_in_cds"
T_RNA_FILE = "t_rna_file"
ASITE_DISP_LENGTH_FILE = "asite_disp_length_file"

# General.
NUM_PROCESSES = "num_processes"
IS_TEST_RUN = "isTestRun"
ALIGNER = "aligner"
