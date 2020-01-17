"""
Test constants.
"""
import os.path
import riboviz

SIMDATA_DIR = os.path.join(riboviz.DATA_DIR, "simdata")
""" Path of data/simdata/ derived from DATA_DIR. """

VIGNETTE_DIR = os.path.join(riboviz.BASE_PATH, "vignette")
""" Path of vignette/ derived from BASE_PATH. """
VIGNETTE_CONFIG = os.path.join(VIGNETTE_DIR, "vignette_config.yaml")
""" Path of vignette/vignette_config.yaml derived from VIGNETTE_DIR. """
VIGNETTE_SAMPLES = ["WT3AT", "WTnone"]
""" Sample names in vignette/vignette_config.yaml. """
VIGNETTE_INDEX_DIR = os.path.join(VIGNETTE_DIR, "index")
""" Path of vignette/index derived from VIGNETTE_DIR. """
VIGNETTE_TMP_DIR = os.path.join(VIGNETTE_DIR, "tmp")
""" Path of vignette/tmp derived from VIGNETTE_DIR. """
VIGNETTE_OUTPUT_DIR = os.path.join(VIGNETTE_DIR, "output")
""" Path of vignette/output derived from VIGNETTE_DIR. """

INDEX_PREFIXES = ["YAL_CDS_w_250", "yeast_rRNA"]
""" Index file prefixes. """
NUM_INDICES = 9
""" Number of index files for each type of index. """

SIMDATA_UMI_CONFIG = os.path.join(VIGNETTE_DIR, "simdata_umi_config.yaml")
"""
Path of vignette/simdata_umi_config.yaml derived from VIGNETTE_DIR.
"""
SIMDATA_UMI_SAMPLE = "umi5_umi3"
""" Sample name in vignette/simdata_umi_config.yaml. """

SIMDATA_MULTIPLEX_CONFIG = os.path.join(
    VIGNETTE_DIR, "simdata_multiplex_config.yaml")
"""
Path of vignette/simdata_multiplex_config.yaml derived from VIGNETTE_DIR.
"""
