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
* ``StopInCDS`` => ``stop_in_feature``
* ``StopInFeature`` => ``stop_in_feature``
* ``codon_pos`` => ``codon_positions_file``
* ``nprocesses`` => ``num_processes``
* ``orf_fasta`` => ``orf_fasta_file``
* ``orf_index`` => ``orf_index_prefix``
* ``ribovizGFF`` => ``is_riboviz_gff``
* ``rRNA_fasta`` => ``rrna_fasta_file``
* ``rRNA_index`` => ``rrna_index_prefix``
* ``stop_in_cds`` => ``stop_in_feature``
* ``t_rna`` => ``t_rna_file``
* ``do_pos_sp_nt_freq`` => ``output_metagene_normalized_profile``

Expected parameters added to the current release are added along
with default values, if they are not already present in the
configuration. These are taken from file
:py:const:`riboviz.params.DEFAULT_CONFIG_YAML`.

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
import riboviz
from riboviz import params


RENAMES = {
    # Names in pre-commit 8da8071, 18 Dec 2019, to current names.
    "Buffer": params.BUFFER,
    "MaxReadLen": params.MAX_READ_LENGTH,
    "MinReadLen": params.MIN_READ_LENGTH,
    "PrimaryID": params.PRIMARY_ID,
    "SecondID": params.SECONDARY_ID,
    "StopInCDS": params.STOP_IN_FEATURE,
    "codon_pos": params.CODON_POSITIONS_FILE,
    "nprocesses": params.NUM_PROCESSES,
    "orf_fasta": params.ORF_FASTA_FILE,
    "orf_index": params.ORF_INDEX_PREFIX,
    "rRNA_fasta": params.RRNA_FASTA_FILE,
    "rRNA_index": params.RRNA_INDEX_PREFIX,
    "ribovizGFF": params.IS_RIBOVIZ_GFF,
    "stop_in_cds": params.STOP_IN_FEATURE,
    "t_rna": params.T_RNA_FILE,
    "do_pos_sp_nt_freq": params.OUTPUT_METAGENE_NORMALIZED_PROFILE
}
"""
Renamed configuration parameters.
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
    configuration. New parameters and default values are taken from
    :py:const:`riboviz.params.DEFAULT_CONFIG_YAML`.

    :param config: Configuration
    :type config: dict
    """
    default_config_file = os.path.join(os.path.dirname(riboviz.__file__),
                                       params.DEFAULT_CONFIG_YAML_FILE)
    with open(default_config_file, "r") as f:
        default_config = yaml.load(f, yaml.SafeLoader)
    # Rename existing parameters.
    for (old_key, new_key) in list(RENAMES.items()):
        if old_key in config:
            value = config[old_key]
            del config[old_key]
            config[new_key] = value
    # Add new parameters.
    for (key, value) in list(default_config.items()):
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
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    else:
        print((yaml.dump(config, sort_keys=False)))
