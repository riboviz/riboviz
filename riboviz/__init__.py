"""
General constants.
"""
import os.path

__version_info__ = ('2', '1')
""" Version """
__version__ = '.'.join(__version_info__)
""" Version """

BASE_PATH = os.path.dirname(os.path.dirname(__file__))
""" Path to parent of :py:mod:`riboviz` module. """
PY_SCRIPTS = os.path.join(BASE_PATH, "riboviz/tools")
""" Path to ``riboviz/tools/`` directory. """
R_SCRIPTS = os.path.join(BASE_PATH, "rscripts")
""" Path to ``rscripts/`` directory. """
