"""
General constants.
"""
import os.path


BASE_PATH = os.path.dirname(os.path.dirname(__file__))
""" Path of parent riboviz module. """
PY_SCRIPTS = os.path.join(BASE_PATH, "riboviz/tools")
""" Path of riboviz/tools/ derived from BASE_PATH. """
R_SCRIPTS = os.path.join(BASE_PATH, "rscripts")
""" Path of rscripts/ derived from BASE_PATH. """
DATA_DIR = os.path.join(BASE_PATH, "data")
""" Path of data/ derived from BASE_PATH. """
