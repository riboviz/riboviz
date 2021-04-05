"""
Test constants.
"""
import os.path
import riboviz

DATA_DIR = os.path.join(riboviz.BASE_PATH, "data")
""" Path to ``data/`` directory. """
SIMDATA_DIR = os.path.join(DATA_DIR, "simdata")
""" Path to ``data/simdata/``. """
VIGNETTE_DIR = os.path.join(riboviz.BASE_PATH, "vignette")
""" Path to ``vignette/``. """
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

SIMDATA_UMI_CONFIG = os.path.join(VIGNETTE_DIR, "simdata_umi_config.yaml")
""" Path to ``vignette/simdata_umi_config.yaml``. """
SIMDATA_UMI_SAMPLE = "umi5_umi3"
""" Sample name in ``vignette/simdata_umi_config.yaml``. """

SIMDATA_MULTIPLEX_CONFIG = os.path.join(
    VIGNETTE_DIR, "simdata_multiplex_config.yaml")
""" Path to ``vignette/simdata_multiplex_config.yaml``. """

NEXTFLOW_WORKFLOW = os.path.join(riboviz.BASE_PATH, "prep_riboviz.nf")
""" Path to ``prep_riboviz.nf``. """
