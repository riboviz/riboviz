"""
Test constants.
"""
import os.path
import riboviz

SIMDATA_DIR = os.path.join(riboviz.DATA_DIR, "simdata")
""" Path to ``data/simdata/``. """
VIGNETTE_DIR = os.path.join(riboviz.BASE_PATH, "vignette")
""" Path to ``vignette/``. """
VIGNETTE_CONFIG = os.path.join(VIGNETTE_DIR, "vignette_config.yaml")
""" Path to ``vignette/vignette_config.yaml``. """
VIGNETTE_SAMPLES = ["WT3AT", "WTnone"]
""" Sample names in ``vignette/vignette_config.yaml``. """
VIGNETTE_INDEX_DIR = os.path.join(VIGNETTE_DIR, "index")
""" Path to ``vignette/index/``. """
VIGNETTE_TMP_DIR = os.path.join(VIGNETTE_DIR, "tmp")
""" Path to ``vignette/tmp/``. """
VIGNETTE_OUTPUT_DIR = os.path.join(VIGNETTE_DIR, "output")
""" Path to ``vignette/output/``. """
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
