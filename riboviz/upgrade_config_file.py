#!/usr/bin/env python
"""
Upgrade YAML configuration file be compatible with current version of
RiboViz.
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
            "StopInCDS": params.STOP_IN_CDS,
            "isTestRun": params.IS_TEST_RUN,
            "ribovizGFF": params.IS_RIBOVIZ_GFF,
            "t_rna": params.T_RNA_FILE,
            "codon_pos": params.CODON_POSITIONS_FILE}
"""
Mapping from parameter names pre-commit 8da8071, 18 Dec 2019, to
current configuration parameter names
"""

UPDATES_10_11 = {
    params.DO_POS_SP_NT_FREQ: True,
    params.FEATURES_FILE:  "data/yeast_features.tsv"
}
"""
Parameters added between release 1.0.0, 9 Oct 2017, 83027ef and 1.1.0,
31 Jan 2019, 340b9b5
"""

UPDATES_11_DEVELOP = {
    params.LOGS_DIR: "vignette/logs",
    params.CMD_FILE: "run_riboviz_vignette.sh",
    params.T_RNA_FILE: "data/yeast_tRNAs.tsv",
    params.CODON_POSITIONS_FILE: "data/yeast_codon_pos_i200.RData",
    params.COUNT_THRESHOLD: 64,
    params.ASITE_DISP_LENGTH_FILE: "data/yeast_standard_asite_disp_length.txt",
    params.COUNT_READS: True
}
"""
Parameters added between release 1.1.0, 31 Jan 2019, 340b9b5 to
pre-commit 8da8071, 18 Dec 2019
"""


def upgrade_config_file(input_file, output_file=None):
    """
    Upgrade YAML configuration file be compatible with current version
    of RiboViz.

    If output_file is None then provided then upgraded content is
    printed to standard out.

    :param input_file: Input YAML file
    :type input_file: str or unicode
    :param output_file: Output YAML file or None
    :type output_file: str or unicode
    :raises AssertionError: if input_file does not exist or is not a
    file
    """
    assert os.path.exists(input_file) and os.path.isfile(input_file),\
        "{} does not exist or is not a file".format(input_file)
    with open(input_file, 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)

    # Upgrade existing keys
    for (old_key, new_key) in list(UPGRADES.items()):
        if old_key in config:
            value = config[old_key]
            del config[old_key]
            config[new_key] = value

    # Parameters added between release 1.0.0, 9 Oct 2017, 83027ef and
    # 1.1.0, 31 Jan 2019, 340b9b5.
    for (key, value) in list(UPDATES_10_11.items()):
        if key not in config:
            config[key] = value

    # Parameters added between release 1.1.0, 31 Jan 2019, 340b9b5 to
    # pre-commit 8da8071, 18 Dec 2019
    for (key, value) in list(UPDATES_11_DEVELOP.items()):
        if key not in config:
            config[key] = value

    # Parameters changed between release 1.1.0, 31 Jan 2019, 340b9b5
    # to pre-commit 8da8071, 18 Dec 2019
    for key in [params.RRNA_INDEX_PREFIX, params.ORF_INDEX_PREFIX]:
        # Index prefixes are now relative to params.DIR_INDEX
        prefix = os.path.split(config[key])[1]
        config[key] = prefix

    # Replace <PATH>/scripts/yeast_features.tsv with
    # <PATH>/data/yeast_features.tsv
    features_file = config[params.FEATURES_FILE]
    features_dir_file = os.path.split(features_file)
    features_super_dir = os.path.split(features_dir_file[0])[0]
    config[params.FEATURES_FILE] = os.path.join(
        features_super_dir, "data", features_dir_file[1])

    if output_file is not None:
        with open(output_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
    else:
        print((yaml.dump(config)))
