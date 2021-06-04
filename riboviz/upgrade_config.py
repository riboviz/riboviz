"""
Upgrade previous versions of the workflow configuration to be
compatible with current version.

Configuration parameters that have been renamed from 1.x are updated:

* ``rRNA_fasta`` => ``rrna_fasta_file``
* ``orf_fasta`` => ``orf_fasta_file``
* ``rRNA_index`` => ``rrna_index_prefix``
* ``orf_index`` => ``orf_index_prefix``
* ``nprocesses`` => ``num_processes``
* ``MinReadLen`` => ``min_read_length``
* ``MaxReadLen`` => ``max_read_length``
* ``Buffer`` => ``buffer``
* ``PrimaryID`` => ``primary_id``
* ``SecondID`` => ``secondary_id``
* ``StopInCDS`` => ``stop_in_feature``
* ``stop_in_cds`` => ``stop_in_feature``
* ``StopInFeature`` => ``stop_in_feature``
* ``ribovizGFF`` => ``is_riboviz_gff``
* ``t_rna`` => ``t_rna_file``
* ``codon_pos`` => ``codon_positions_file``

Expected parameters added between release 1.0.0, 9 Oct 2017, 83027ef
and 1.1.0, 31 Jan 2019, 340b9b5, are added along with default values,
if they are not already present in the configuration:

* ``do_pos_sp_nt_freq: true``
* ``features_file: data/yeast_features.tsv``

Expected parameters added between release 1.1.0 and 2.0 are added
along with default values, if they are not already present in the
configuration:

* ``asite_disp_length_file: data/yeast_standard_asite_disp_length.txt``
* ``cmd_file: run_riboviz_vignette.sh``
* ``codon_positions_file: data/yeast_codon_pos_i200.RData``
* ``count_reads: true``
* ``count_threshold: 64``
* ``dedup_stats: false``
* ``dedup_umis: false``
* ``dir_logs: vignette/logs``
* ``extract_umis: false``
* ``group_umis: false``
* ``multiplex_fq_files: null``
* ``publish_index_tmp: false``
* ``sample_sheet: null``
* ``trim_5p_mismatches: true``
* ``t_rna_file: data/yeast_tRNAs.tsv``
* ``umi_regexp: null``

Expected parameters added between release 2.0 and the current
release are added along with default values, if they are not
already present in the configuration:

* ``feature: CDS ``
* ``job_email: null``
* ``job_email_events: beas``
* ``job_memory: 8G``
* ``job_name: riboviz``
* ``job_num_cpus: 4``
* ``job_parallel_env: mpi``
* ``job_runtime: '48:00:00'``
* ``nextflow_dag_file: nextflow-dag.html``
* ``nextflow_report_file: nextflow-report.html``
* ``nextflow_timeline_file: nextflow-timeline.html``
* ``nextflow_trace_file: nextflow-trace.tsv``
* ``nextflow_work_dir: work``
* ``stop_in_feature: false``
* ``validate_only: false``
* ``run_static_html: true``

The values of parameters ``rrna_index_prefix`` and
``orf_index_prefix`` are updated to be file names only, as, these are
now assumed to be relative to ``<dir_index>``. For example the
configuration parameters::

    rRNA_index: vignette/index/yeast_rRNA
    orf_index: vignette/index/YAL_CDS_w_250

are updated to::

    rRNA_index_prefix: yeast_rRNA
    orf_index_prefix: YAL_CDS_w_250
"""
import os
import os.path
import yaml
from riboviz import params

UPGRADES = {"rRNA_fasta": params.RRNA_FASTA_FILE,
            "orf_fasta": params.ORF_FASTA_FILE,
            "rRNA_index": params.RRNA_INDEX_PREFIX,
            "orf_index": params.ORF_INDEX_PREFIX,
            "nprocesses": params.NUM_PROCESSES,
            "MinReadLen": params.MIN_READ_LENGTH,
            "MaxReadLen": params.MAX_READ_LENGTH,
            "Buffer": params.BUFFER,
            "PrimaryID": params.PRIMARY_ID,
            "SecondID": params.SECONDARY_ID,
            "StopInCDS": params.STOP_IN_FEATURE,
            params.STOP_IN_CDS: params.STOP_IN_FEATURE,
            "ribovizGFF": params.IS_RIBOVIZ_GFF,
            "t_rna": params.T_RNA_FILE,
            "codon_pos": params.CODON_POSITIONS_FILE}
"""
Map from configuration parameter names pre-commit 8da8071, 18 Dec
2019, to current configuration parameter names.
"""

UPDATES_10_11 = {
    params.DO_POS_SP_NT_FREQ: True,
    params.FEATURES_FILE:  "data/yeast_features.tsv"
}
"""
Map from configuration parameters to default values for parameters
added between release 1.0.0 and 1.1.0.
"""

UPDATES_11_20 = {
    params.ASITE_DISP_LENGTH_FILE: "data/yeast_standard_asite_disp_length.txt",
    params.CMD_FILE: "run_riboviz_vignette.sh",
    params.CODON_POSITIONS_FILE: "data/yeast_codon_pos_i200.RData",
    params.COUNT_READS: True,
    params.COUNT_THRESHOLD: 64,
    params.DEDUP_STATS: False,
    params.DEDUP_UMIS: False,
    params.EXTRACT_UMIS: False,
    params.GROUP_UMIS: False,
    params.LOGS_DIR: "vignette/logs",
    params.MULTIPLEX_FQ_FILES: None,
    params.PUBLISH_INDEX_TMP: False,
    params.SAMPLE_SHEET: None,
    params.TRIM_5P_MISMATCHES: True,
    params.T_RNA_FILE: "data/yeast_tRNAs.tsv",
    params.UMI_REGEXP: None
}
"""
Map from configuration parameters to default values for parameters
added between release 1.1.0 and 2.0.
"""

UPDATES_20_CURRENT = {
    params.FEATURE: "CDS",
    params.RUN_STATIC_HTML: True,
    params.STOP_IN_FEATURE: False
}

"""
Map from configuration parameters to default values for parameters
added between release 2.0 and the current release.
"""


def upgrade_config(config):
    """
    Upgrade workflow configuration to be compatible with current
    configuration.

    :param config: Configuration
    :type config: dict
    """
    # Upgrade existing keys
    for (old_key, new_key) in list(UPGRADES.items()):
        if old_key in config:
            value = config[old_key]
            del config[old_key]
            config[new_key] = value

    for (key, value) in list(UPDATES_10_11.items()):
        if key not in config:
            config[key] = value

    for (key, value) in list(UPDATES_11_20.items()):
        if key not in config:
            config[key] = value

    # Remove params.NEXTFLOW_RESUME as it is a command-line only
    # configuration parameter.
    job_config = params.DEFAULT_JOB_CONFIG.copy()
    del job_config[params.NEXTFLOW_RESUME]
    for (key, value) in list(job_config.items()):
        if key not in config:
            config[key] = value

    for (key, value) in list(UPDATES_20_CURRENT.items()):
        if key not in config:
            config[key] = value

    # Parameters changed between release 1.1.0, 31 Jan 2019, 340b9b5
    # to pre-commit 8da8071, 18 Dec 2019
    for key in [params.RRNA_INDEX_PREFIX, params.ORF_INDEX_PREFIX]:
        # Index prefixes are now relative to params.DIR_INDEX
        prefix = os.path.split(config[key])[1]
        config[key] = prefix


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
