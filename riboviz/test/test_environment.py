"""
:py:mod:`riboviz.environment` tests.
"""
import pytest
from riboviz import environment
from riboviz import params


@pytest.mark.parametrize(
    "envs_values",
    [{},
     {params.ENV_RIBOVIZ_SAMPLES: "samples"},
     {params.ENV_RIBOVIZ_ORGANISMS: "organisms"},
     {params.ENV_RIBOVIZ_DATA: "data"},
     {params.ENV_RIBOVIZ_SAMPLES: "samples",
      params.ENV_RIBOVIZ_ORGANISMS: "organisms",
      params.ENV_RIBOVIZ_DATA: "data"}],
    ids=["none", "samples", "organisms", "data", "all"])
def test_get_environment_vars(envs_values, monkeypatch):
    """
    Test :py:func:`riboviz.environment.get_environment_vars`.
    :param envs_values: Dictionary of environment variables to values
    :type envs_values: dict
    :param monkeypatch: Fixture to help modify OS environment variables \
    for purposes of test (pytest built-in fixture)
    :type monkeypatch: _pytest.monkeypatch.MonkeyPatch
    """
    for e, v in envs_values.items():
        monkeypatch.setenv(e, v)
    actual_envs_values = environment.get_environment_vars()
    for e, v in envs_values.items():
        assert e in actual_envs_values.keys(), \
            "Missing value for environment variable {}".format(e)
        assert actual_envs_values[e] == v, \
            "Unexpected nexpected value for environment variable {}".format(e)


@pytest.mark.parametrize("env", params.ENV_DIRS)
@pytest.mark.parametrize("param", params.ENV_PARAMS)
def test_apply_env_to_config(env, param, monkeypatch):
    """
    Test :py:func:`riboviz.environment.apply_env_to_config` with
    a combination of environment variable and RiboViz parameter
    that supports the use of environment variable tokens
    :param env: Environment variable
    :type envs_values: str or unicode
    :param param: Parameter
    :type param: str or unicode
    :param monkeypatch: Fixture to help modify OS environment variables \
    for purposes of test (pytest built-in fixture)
    :type monkeypatch: _pytest.monkeypatch.MonkeyPatch
    """
    monkeypatch.setenv(env, "/home/test")
    # Set params.ADAPTERS too for a configuration parameter that
    # won't be updated.
    config = {
        params.ADAPTERS: "CTGTAGGCACC",
        param: "{}/subdirectory/somefile.txt".format(
            environment.ENV_TOKEN_FORMAT.format(env))
    }
    environment.apply_env_to_config(config)
    expected_config = {
        params.ADAPTERS: "CTGTAGGCACC",
        param: "/home/test/subdirectory/somefile.txt"
    }
    assert config == expected_config, "Unexpected configuration"


@pytest.mark.parametrize("env", params.ENV_DIRS)
@pytest.mark.parametrize("param", params.ENV_PARAMS)
def test_apply_default_env_to_config(env, param):
    """
    Test :py:func:`riboviz.environment.apply_env_to_config` with
    a combination of environment variable and RiboViz parameter
    that supports the use of environment variable tokens. In this
    test the environment variable token is not used so the
    default value of "." is expected to be substituted.
    :param env: Environment variable
    :type envs_values: str or unicode
    :param param: Parameter
    :type param: str or unicode
    :param monkeypatch: Fixture to help modify OS environment variables \
    for purposes of test (pytest built-in fixture)
    :type monkeypatch: _pytest.monkeypatch.MonkeyPatch
    """
    config = {
        params.ADAPTERS: "CTGTAGGCACC",
        param: "{}/subdirectory/somefile.txt".format(
            environment.ENV_TOKEN_FORMAT.format(env))
    }
    environment.apply_env_to_config(config)
    expected_config = {
        params.ADAPTERS: "CTGTAGGCACC",
        param: "{}/subdirectory/somefile.txt".format(
            environment.DEFAULT_ENV_DIR)
    }
    assert config == expected_config, "Unexpected configuration"
