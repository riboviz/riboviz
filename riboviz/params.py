"""
RiboViz configuration parameters.
"""

# Input and output directories.
CMD_FILE = "cmd_file"
INPUT_DIR ="dir_in"
INDEX_DIR ="dir_index"
LOGS_DIR ="dir_logs"
OUTPUT_DIR ="dir_out"
TMP_DIR ="dir_tmp"

# Input files.
FEATURES_FILE = "features_file"
ORF_FASTA_FILE = "orf_fasta"
ORF_GFF_FILE = "orf_gff_file"
R_RNA_FASTA_FILE = "rRNA_fasta"

# Indexing.
BUILD_INDICES = "build_indices"
ORF_INDEX_PREFIX = "orf_index"
R_RNA_INDEX_PREFIX = "rRNA_index"

# Adapter trimming.
ADAPTERS = "adapters"

# Sample files.
FQ_FILES = "fq_files"

# Demultiplexing and deduplication.
SAMPLE_SHEET_FILE = "sample_sheet"
DEDUP_UMIS = "dedup_umis"
GROUP_UMIS = "group_umis"
EXTRACT_UMIS = "extract_umis"
MULTIPLEX_FQ_FILES = "multiplex_fq_files"
UMI_REGEXP = "umi_regexp"

# Bedgraphs.
MAKE_BEDGRAPH = "make_bedgraph"

# Bam to H5 and statistics generation.
BUFFER = "Buffer"
CODON_POS = "codon_pos"
COUNT_THRESHOLD = "count_threshold"
DATASET = "dataset"
DO_POS_SP_NT_FREQ = "do_pos_sp_nt_freq"
MAX_READ_LEN = "MaxReadLen"
MIN_READ_LEN = "MinReadLen"
PRIMARY_ID = "PrimaryID"
RIBOVIZ_GFF = "ribovizGFF"
RPF = "rpf"
SECOND_ID = "SecondID"
STOP_IN_CDS = "StopInCDS"
T_RNA = "t_rna"

# General.
NPROCESSES = "nprocesses"
IS_TEST_RUN = "isTestRun"
ALIGNER = "aligner"
