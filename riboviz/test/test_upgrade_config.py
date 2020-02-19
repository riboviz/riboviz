"""
:py:mod:`riboviz.upgrade_config` tests.
"""
import os
import tempfile
import pytest
import yaml
from riboviz import upgrade_config
import riboviz.test
from riboviz.test import config


@pytest.fixture(scope="function")
def tmp_file():
    """
    Create a temporary file with a ``yaml`` suffix.

    :return: path to temporary file
    :rtype: str or unicode
    """
    _, tmp_file = tempfile.mkstemp(prefix="tmp", suffix=".yaml")
    yield tmp_file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


def upgrade_and_validate(config_file,
                         upgraded_config_file,
                         expected_config_file):
    """
    Test :py:func:`riboviz.upgrade_config.upgrade_config_file` on a
    configuration file and compare to expected contents.

    :param config_file: configuration file
    :type config_file: str or unicode
    :param upgraded_config_file: upgraded configuration file
    :type upgraded_config_file: str or unicode
    :param expected_config_file: expected configuration file
    :type expected_config_file: str or unicode
    :raises AssertionError: if contents differ
    """
    config_path = os.path.join(os.path.dirname(config.__file__),
                               config_file)
    upgrade_config.upgrade_config_file(config_path, upgraded_config_file)
    with open(expected_config_file, "r") as f:
        expected_config = yaml.load(f, yaml.SafeLoader)
    with open(upgraded_config_file, "r") as f:
        actual_config = yaml.load(f, yaml.SafeLoader)
    assert expected_config == actual_config,\
        "Upgraded {} differs from {}".format(config_file,
                                             expected_config_file)


@pytest.mark.parametrize("config_file",
                         ["vignette_config_1.0.0.yaml",
                          "vignette_config_1.1.0.yaml",
                          "vignette_config_current.yaml"])
def test_upgrade_vignette_config(config_file, tmp_file):
    """
    Test :py:func:`riboviz.upgrade_config.upgrade_config_file` on
    previous versions of the vignette configuration file
    (in :py:mod:`riboviz.test.config`) and compare to the current
    version.

    :param config_file: configuration file
    :type config_file: str or unicode
    """
    upgrade_and_validate(config_file, tmp_file,
                         riboviz.test.VIGNETTE_CONFIG)


@pytest.mark.parametrize("config_file",
                         ["simdata_umi_config_current.yaml"])
def test_upgrade_simdata_umi_config(config_file, tmp_file):
    """
    Test :py:func:`riboviz.upgrade_config.upgrade_config_file` on
    previous versions of the simulated UMI data configuration file
    (in :py:mod:`riboviz.test.config`) and compare to the current
    version.

    :param config_file: configuration file
    :type config_file: str or unicode
    :param tmp_file: upgraded configuration file
    :type tmp_file: str or unicode
    """
    upgrade_and_validate(config_file, tmp_file,
                         riboviz.test.SIMDATA_UMI_CONFIG)


@pytest.mark.parametrize("config_file",
                         ["simdata_multiplex_config_current.yaml"])
def test_upgrade_simdata_multiplex_config(config_file, tmp_file):
    """
    Test :py:func:`riboviz.upgrade_config.upgrade_config_file` on
    previous versions of the simulated multiplexed data configuration
    file (in :py:mod:`riboviz.test.config`) and compare to the current
    version.

    :param config_file: configuration file
    :type config_file: str or unicode
    :param tmp_file: upgraded configuration file
    :type tmp_file: str or unicode
    """
    upgrade_and_validate(config_file, tmp_file,
                         riboviz.test.SIMDATA_MULTIPLEX_CONFIG)
