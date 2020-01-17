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
SIMDATA_UMI_CONFIG = os.path.join(VIGNETTE_DIR, "simdata_umi_config.yaml")
"""
Path of vignette/simdata_umi_config.yaml derived from VIGNETTE_DIR.
"""
SIMDATA_MULTIPLEX_CONFIG = os.path.join(
    VIGNETTE_DIR, "simdata_multiplex_config.yaml")
"""
Path of vignette/simdata_multiplex_config.yaml derived from VIGNETTE_DIR.
"""
