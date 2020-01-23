"""
riboviz.upgrade_config_file test suite.
"""
import os
import tempfile
import pytest
import yaml
import riboviz.test
from riboviz.test import config
from riboviz.upgrade_config_file import upgrade_config_file


@pytest.fixture(scope="function")
def upgraded_config_file():
    """
    Create a temporary file to write upgraded configuration into.

    :return: path to configuration file
    :rtype: str or unicode
    """
    _, config_file = tempfile.mkstemp(prefix="tmp", suffix=".yaml")
    yield config_file
    if os.path.exists(config_file):
        os.remove(config_file)


def upgrade_and_validate(config_file,
                         upgraded_config_file,
                         expected_config_file):
    """
    Run upgrade_config_file on a configuration file and compare to
    expected contents.

    :param config_file: configuration file to upgrade
    :type config_file: str or unicode
    :param upgraded_config_file: path to write upgraded configuration
    file into (pytest fixture)
    :type upgraded_config_file: str or unicode
    :param expected_config_file: configuration file with expectec
    contents
    :type expected_config_file: str or unicode
    :raises AssertionError: if contents differ
    """
    config_path = os.path.join(os.path.dirname(config.__file__),
                               config_file)
    upgrade_config_file(config_path, upgraded_config_file)
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
                          "vignette_config_a5bccc6_20200123.yaml"])
def test_upgrade_vignette_config(config_file, upgraded_config_file):
    """
    Run upgrade_config_file on previous versions of the vignette
    configuration file and compare to the current version.

    :param config_file: configuration file to upgrade
    :type config_file: str or unicode
    :param upgraded_config_file: path to write upgraded configuration
    file into (pytest fixture)
    :type upgraded_config_file: str or unicode
    """
    upgrade_and_validate(config_file,
                         upgraded_config_file,
                         riboviz.test.VIGNETTE_CONFIG)


@pytest.mark.parametrize("config_file",
                         ["simdata_umi_config_a5bccc6_20200123.yaml"])
def test_upgrade_simdata_umi_config(config_file, upgraded_config_file):
    """
    Run upgrade_config_file on previous versions of the simulated UMI
    data configuration file and compare to the current version.

    :param config_file: configuration file to upgrade
    :type config_file: str or unicode
    :param upgraded_config_file: path to write upgraded configuration
    file into (pytest fixture)
    :type upgraded_config_file: str or unicode
    """
    upgrade_and_validate(config_file,
                         upgraded_config_file,
                         riboviz.test.SIMDATA_UMI_CONFIG)


@pytest.mark.parametrize("config_file",
                         ["simdata_multiplex_config_a5bccc6_20200123.yaml"])
def test_upgrade_simdata_multiplex_config(config_file, upgraded_config_file):
    """
    Run upgrade_config_file on previous versions of the simulated
    multiplexed data configuration file and compare to the current
    version.

    :param config_file: configuration file to upgrade
    :type config_file: str or unicode
    :param upgraded_config_file: path to write upgraded configuration
    file into (pytest fixture)
    :type upgraded_config_file: str or unicode
    """
    upgrade_and_validate(config_file,
                         upgraded_config_file,
                         riboviz.test.SIMDATA_MULTIPLEX_CONFIG)
