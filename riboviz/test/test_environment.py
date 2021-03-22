"""
:py:mod:`riboviz.environment` tests.
"""
import os
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
def test_get_environment_vars(envs_values):
    """
    Test :py:func:`riboviz.environment.get_environment_vars`.
    :param envs_values: Dictionary of environment variables to values
    :type envs_values: dict
    """
    os.environ.update(envs_values)
    actual_envs_values = environment.get_environment_vars()
    for e, v in envs_values.items():
        assert e in actual_envs_values.keys(), \
            "Missing value for environment variable {}".format(e)
        assert actual_envs_values[e] == v, \
            "Unexpected nexpected value for environment variable {}".format(e)
