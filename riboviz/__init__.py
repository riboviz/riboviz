"""
General constants.
"""
import os.path
import riboviz


BAM_TO_H5_R = "bam_to_h5.R"
""" Name of bam_to_h5.R script """
GENERATE_STATS_FIGS_R = "generate_stats_figs.R"
""" Name of generate_stats_figs.R script """
COLLATE_TPMS_R = "collate_tpms.R"
""" NAme of collate_tpms.R script """

BASE_PATH = os.path.dirname(os.path.dirname(riboviz.__file__))
""" Path of parent riboviz module. """
PY_SCRIPTS = os.path.join(BASE_PATH, "riboviz/tools")
""" Path of riboviz/tools/ derived from BASE_PATH. """
R_SCRIPTS = os.path.join(BASE_PATH, "rscripts")
""" Path of rscripts/ derived from BASE_PATH. """
DATA_DIR = os.path.join(BASE_PATH, "data")
""" Path of data/ derived from BASE_PATH. """
SIMDATA_DIR = os.path.join(DATA_DIR, "simdata")
""" Path of data/simdata/ derived from DATA_DIR. """

VIGNETTE_DIR = os.path.join(BASE_PATH, "vignette")
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
