"""
General test-related constants.
"""
import os.path
import riboviz


BASE_PATH = os.path.dirname(os.path.dirname(riboviz.__file__))
""" Path of parent riboviz module. """

VIGNETTE_DIR = os.path.join(BASE_PATH, "vignette")
""" Path of vignette/ derived from BASE_PATH. """

PY_SCRIPTS = os.path.join(BASE_PATH, "riboviz/tools")
""" Path of riboviz/tools/ derived from BASE_PATH. """

R_SCRIPTS = os.path.join(BASE_PATH, "rscripts")
""" Path of rscripts/ derived from BASE_PATH. """

DATA_DIR = os.path.join(BASE_PATH, "data")
""" Path of data/ derived from BASE_PATH. """

EXAMPLE_DATA_DIR = os.path.join(DATA_DIR, "example")
""" Path of data/example/ derived from DATA_DIR. """

VIGNETTE_CONFIG = os.path.join(VIGNETTE_DIR, "vignette_config.yaml")
""" Path of vignette/vignette_config.yaml derived from VIGNETTE_DIR. """

EXAMPLE_CONFIG = os.path.join(VIGNETTE_DIR, "example_config.yaml")
"""
Path of vignette/example_config.yaml derived from VIGNETTE_DIR.
"""
