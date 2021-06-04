"""
Upgrade previous versions of the workflow configuration to be
compatible with current version.

Configuration parameters that have been renamed are updated (their
existing values are preserved):

* ``Buffer`` => ``buffer``
* ``MaxReadLen`` => ``max_read_length``
* ``MinReadLen`` => ``min_read_length``
* ``PrimaryID`` => ``primary_id``
* ``SecondID`` => ``secondary_id``
* ``StopInCDS`` => ``stop_in_cds``
* ``codon_pos`` => ``codon_positions_file``
* ``nprocesses`` => ``num_processes``
* ``orf_fasta`` => ``orf_fasta_file``
* ``orf_index`` => ``orf_index_prefix``
* ``rRNA_fasta`` => ``rrna_fasta_file``
* ``rRNA_index`` => ``rrna_index_prefix``
* ``ribovizGFF`` => ``is_riboviz_gff``
* ``t_rna`` => ``t_rna_file``

Expected parameters added to the current release are added along
with default values, if they are not already present in the
configuration:

* ``asite_disp_length_file: data/yeast_standard_asite_disp_length.txt``
* ``codon_positions_file: data/yeast_codon_pos_i200.RData``
* ``count_reads: true``
* ``count_threshold: 64``
* ``dedup_stats: false``
* ``dedup_umis: false``
* ``do_pos_sp_nt_freq: true``
* ``extract_umis: false``
* ``features_file: data/yeast_features.tsv``
* ``group_umis: false``
* ``job_email: null``
* ``job_email_events: beas``
* ``job_memory: 8G``
* ``job_name: riboviz``
* ``job_num_cpus: 4``
* ``job_parallel_env: mpi``
* ``job_runtime: '48:00:00'``
* ``multiplex_fq_files: null``
* ``nextflow_dag_file: nextflow-dag.html``
* ``nextflow_report_file: nextflow-report.html``
* ``nextflow_timeline_file: nextflow-timeline.html``
* ``nextflow_trace_file: nextflow-trace.tsv``
* ``nextflow_work_dir: work``
* ``publish_index_tmp: false``
* ``sample_sheet: null``
* ``samsort_memory: null``
* ``t_rna_file: data/yeast_tRNAs.tsv``
* ``trim_5p_mismatches: true``
* ``umi_regexp: null``
* ``validate_only: false``
* ``run_static_html: true``

The values of parameters ``rrna_index_prefix`` and
``orf_index_prefix`` are updated to be file names only, as, these are
now assumed to be relative to ``<dir_index>``. For example the
configuration parameters::

    rRNA_index: vignette/index/yeast_rRNA
    orf_index: vignette/index/YAL_CDS_w_250

are updated to::

    rrna_index_prefix: yeast_rRNA
    orf_index_prefix: YAL_CDS_w_250

Configuration parameters that are now unused are removed:

* ``aligner``
* ``isTestRun``
* ``is_test_run``
* ``cmd_file``
* ``dir_logs``
"""
import os
import os.path
import yaml
from riboviz import params


RENAMES = {
    # Names in pre-commit 8da8071, 18 Dec 2019, to current names.
    "Buffer": params.BUFFER,
    "MaxReadLen": params.MAX_READ_LENGTH,
    "MinReadLen": params.MIN_READ_LENGTH,
    "PrimaryID": params.PRIMARY_ID,
    "SecondID": params.SECONDARY_ID,
    "StopInCDS": params.STOP_IN_CDS,
    "codon_pos": params.CODON_POSITIONS_FILE,
    "nprocesses": params.NUM_PROCESSES,
    "orf_fasta": params.ORF_FASTA_FILE,
    "orf_index": params.ORF_INDEX_PREFIX,
    "rRNA_fasta": params.RRNA_FASTA_FILE,
    "rRNA_index": params.RRNA_INDEX_PREFIX,
    "ribovizGFF": params.IS_RIBOVIZ_GFF,
    "t_rna": params.T_RNA_FILE
}
"""
Renamed configuration parameters.
"""

UPDATES = {
    params.ASITE_DISP_LENGTH_FILE: "data/yeast_standard_asite_disp_length.txt",
    params.DO_POS_SP_NT_FREQ: True,
    params.CODON_POSITIONS_FILE: "data/yeast_codon_pos_i200.RData",
    params.COUNT_READS: True,
    params.COUNT_THRESHOLD: 64,
    params.DEDUP_STATS: False,
    params.DEDUP_UMIS: False,
    params.EXTRACT_UMIS: False,
    params.FEATURES_FILE:  "data/yeast_features.tsv",
    params.FQ_FILES: None,
    params.GROUP_UMIS: False,
    params.MULTIPLEX_FQ_FILES: None,
    params.PUBLISH_INDEX_TMP: False,
    params.RUN_STATIC_HTML: True,
    params.SAMPLE_SHEET: None,
    params.SAMSORT_MEMORY: None,
    params.TRIM_5P_MISMATCHES: True,
    params.T_RNA_FILE: "data/yeast_tRNAs.tsv",
    params.UMI_REGEXP: None
}
"""
Map from configuration parameters to default values for parameters.
"""

UNUSED = [
    "aligner",
    "isTestRun",
    "is_test_run",
    "cmd_file",
    "dir_logs"
]
"""
Unused configuration parameters for removal.
"""


def upgrade_config(config):
    """
    Upgrade workflow configuration to be compatible with current
    configuration.

    :param config: Configuration
    :type config: dict
    """
    # Rename existing parameters.
    for (old_key, new_key) in list(RENAMES.items()):
        if old_key in config:
            value = config[old_key]
            del config[old_key]
            config[new_key] = value

    # Add new parameters.
    for (key, value) in list(UPDATES.items()):
        if key not in config:
            config[key] = value
    # Remove params.NEXTFLOW_RESUME as it is a command-line only
    # configuration parameter.
    job_config = params.DEFAULT_JOB_CONFIG.copy()
    del job_config[params.NEXTFLOW_RESUME]
    for (key, value) in list(job_config.items()):
        if key not in config:
            config[key] = value

    # Index prefixes are now relative to params.DIR_INDEX
    for key in [params.RRNA_INDEX_PREFIX, params.ORF_INDEX_PREFIX]:
        prefix = os.path.split(config[key])[1]
        config[key] = prefix

    # Removed unused parameters.
    for key in UNUSED:
        if key in config:
            del config[key]


def upgrade_config_file(input_file, output_file=None):
    """
    Upgrade workflow configuration file to be compatible with current
    configuration.

    If ``output_file`` is ``None`` then the upgraded configuration is
    printed to standard output.

    :param input_file: Input file
    :type input_file: str or unicode
    :param output_file: Output file or None
    :type output_file: str or unicode
    :raises AssertionError: If ``input_file`` does not exist or is \
    not a file
    """
    assert os.path.exists(input_file) and os.path.isfile(input_file),\
        "{} does not exist or is not a file".format(input_file)
    with open(input_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    upgrade_config(config)
    if output_file is not None:
        with open(output_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
    else:
        print((yaml.dump(config)))
