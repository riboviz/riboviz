"""
Environment-related functions.
"""
import os
from riboviz import params
from riboviz import utils


def get_environment_vars():
    """
    Get map from RiboViz environment variables, to their values,
    if defined. The environment variables are:

    * :py:const:`riboviz.params.ENV_RIBOVIZ_SAMPLES`
    * :py:const:`riboviz.params.ENV_RIBOVIZ_ORGANISMS`
    * :py:const:`riboviz.params.ENV_RIBOVIZ_DATA`

    :return: Map from environment variables to values
    :rtype: dict
    """
    return {env: os.environ[env] for env in
            [params.ENV_RIBOVIZ_SAMPLES,
             params.ENV_RIBOVIZ_ORGANISMS,
             params.ENV_RIBOVIZ_DATA] if env in os.environ}


def apply_env_to_config(config):
    """
    Replace environment variable tokens with environment variables
    in configuration parameter values that support environment
    variables. See :py:const:`riboviz.params.ENV_PARAMS` for the
    relevant parameters. All other parameters are left unchanged.

    :param config: Configuration
    :type config: dict
    """
    env_vars = get_environment_vars()
    env_tokens = {"${{{0}}}".format(env): value
                  for env, value in env_vars.items()}
    for param in params.ENV_PARAMS:
        if param in config and config[param]:
            config[param] = utils.replace_tokens(
                config[param], env_tokens)
