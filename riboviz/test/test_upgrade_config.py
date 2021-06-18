"""
:py:mod:`riboviz.upgrade_config` tests.
"""
import os
import tempfile
import pytest
import yaml
from riboviz import test
from riboviz import upgrade_config
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
    expected_config_path = os.path.join(
        os.path.dirname(config.__file__),
        expected_config_file)
    with open(expected_config_path, "r") as f:
        expected_config = yaml.load(f, yaml.SafeLoader)
    with open(upgraded_config_file, "r") as f:
        actual_config = yaml.load(f, yaml.SafeLoader)
    assert expected_config == actual_config,\
        "Upgraded {} differs from {}".format(config_file,
                                             expected_config_file)


@pytest.mark.parametrize("config_file",
                         ["config_1.0.0.yaml",
                          "config_1.1.0.yaml",
                          "config_2.0.0.yaml",
                          "config_current.yaml"])
def test_upgrade_config_file(config_file, tmp_file):
    """
    Test :py:func:`riboviz.upgrade_config.upgrade_config_file` on
    previous versions of the configuration file
    (in :py:mod:`riboviz.test.config`) and compare to the current
    version.

    :param config_file: configuration file
    :type config_file: str or unicode
    """
    upgrade_and_validate(config_file, tmp_file,
                         "config_current.yaml")


@pytest.mark.parametrize("config_file",
                         [test.VIGNETTE_CONFIG,
                          test.SIMDATA_UMI_CONFIG,
                          test.SIMDATA_MULTIPLEX_CONFIG])
def test_upgrade_config_file_identity(config_file, tmp_file):
    """
    Test :py:func:`riboviz.upgrade_config.upgrade_config_file` on
    current versions of configuration files and ensure that the
    result is the same i.e. the configurations are already up
    to date and no parameters are changed during the upgrade.

    :param config_file: configuration file
    :type config_file: str or unicode
    """
    upgrade_and_validate(config_file, tmp_file, config_file)
