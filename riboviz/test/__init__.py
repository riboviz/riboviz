"""
Test constants and functions.
"""
import os
import os.path
import riboviz

DATA_DIR = os.path.join(riboviz.BASE_PATH, "data")
""" Path to ``data/`` directory. """
SIMDATA_DIR = os.path.join(DATA_DIR, "simdata")
""" Path to ``data/simdata/``. """
VIGNETTE_DIR = os.path.join(riboviz.BASE_PATH, "vignette")
""" Path to ``vignette/``. """
VIGNETTE_INPUT_DIR = os.path.join(VIGNETTE_DIR, "input")
""" Path to ``vignette/input/``. """
VIGNETTE_CONFIG = os.path.join(VIGNETTE_DIR, "vignette_config.yaml")
""" Path to ``vignette/vignette_config.yaml``. """
VIGNETTE_SAMPLES = ["WT3AT", "WTnone"]
""" Sample names in ``vignette/vignette_config.yaml``. """
VIGNETTE_MISSING_SAMPLE = "NotHere"
"""
Name of missing sample in ``vignette/vignette_config.yaml`` for tests.
"""
INDEX_PREFIXES = ["YAL_CDS_w_250", "yeast_rRNA"]
""" Index file prefixes. """
NUM_INDICES = 9
""" Number of index files for each type of index. """
ORGANISM_FILES = [
    os.path.join(VIGNETTE_DIR, "input", f) for f in [
        "yeast_rRNA_R64-1-1.fa",
        "yeast_YAL_CDS_w_250utrs.fa",
        "yeast_YAL_CDS_w_250utrs.gff3"
    ]
]
""" Organism files in :py:const:`VIGNETTE_INPUT_DIR`. """
DATA_FILES = [
    os.path.join(DATA_DIR, f) for f in [
        "yeast_codon_pos_i200.RData",
        "yeast_features.tsv",
        "yeast_standard_asite_disp_length.txt",
        "yeast_tRNAs.tsv"
    ]
]
""" Data files in :py:const:`DATA_DIR`. """
VIGNETTE_INPUT_FILES = [
    os.path.join(VIGNETTE_INPUT_DIR, f) for f in [
        "SRR1042855_s1mi.fastq.gz",
        "SRR1042864_s1mi.fastq.gz"
        ]
]
""" Vignette input files in :py:const:`VIGNETTE_INPUT_DIR`. """
SIMDATA_INPUT_FILES = [
    os.path.join(SIMDATA_DIR, f) for f in [
        "umi5_umi3_umi_adaptor.fastq",
        "multiplex_umi_barcode_adaptor.fastq",
        "multiplex_barcodes.tsv"
    ]
]
""" Simulated data input files in :py:const:`SIMDATA_DIR`. """

SIMDATA_UMI_CONFIG = os.path.join(VIGNETTE_DIR, "simdata_umi_config.yaml")
""" Path to ``vignette/simdata_umi_config.yaml``. """
SIMDATA_UMI_SAMPLE = "umi5_umi3"
""" Sample name in ``vignette/simdata_umi_config.yaml``. """

SIMDATA_MULTIPLEX_CONFIG = os.path.join(
    VIGNETTE_DIR, "simdata_multiplex_config.yaml")
""" Path to ``vignette/simdata_multiplex_config.yaml``. """

NEXTFLOW_WORKFLOW = os.path.join(riboviz.BASE_PATH, "prep_riboviz.nf")
""" Path to ``prep_riboviz.nf``. """


def customise_path(prefix, path):
    """
    Customise a file path. Get the basename of the path, and create a
    new path, formed from that and ``prefix``.

    If the original path ends with ``/`` then this is first removed to
    allow the correct basepath to be extracted.

    Examples, assuming prefix is ``${RIBOVIZ_SAMPLES}``:

    * ``data/simdata/`` => ``${RIBOVIZ_SAMPLES}/simdata``
    * ``data/simdata`` => ``${RIBOVIZ_SAMPLES}/simdata``

    :param prefix: Prefix
    :type prefix: str or unicode
    :param path: Pathe
    :type path: str or unicode
    :return: Updated pathe
    :rtype: str or unicodee
    """
    nu_path = path
    if path:
        if path.endswith("/"):
            nu_path = path[:-1]
        nu_path = os.path.join(prefix, os.path.basename(nu_path))
    return nu_path


def symlink_files(directory, files):
    """
    Create symbolic links in a directory to each of a list of files.

    :param directory: Directory in which to create symbolic links
    :type directory: py._path.local.LocalPath or str or unicode
    :param files: Files
    :type files: list(str or unicode)
    """
    for f in files:
        os.symlink(f, os.path.join(directory, os.path.basename(f)))
