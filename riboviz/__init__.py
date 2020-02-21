"""
General constants.
"""
import os.path

BASE_PATH = os.path.dirname(os.path.dirname(__file__))
""" Path to parent of :py:mod:`riboviz` module. """
PY_SCRIPTS = os.path.join(BASE_PATH, "riboviz/tools")
""" Path to ``riboviz/tools/`` directory. """
R_SCRIPTS = os.path.join(BASE_PATH, "rscripts")
""" Path to ``rscripts/`` directory. """
DATA_DIR = os.path.join(BASE_PATH, "data")
""" Path to ``data/`` directory. """
