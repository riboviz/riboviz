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
