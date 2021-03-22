"""
Environment-related functions.
"""
import os
from riboviz import params
from riboviz import utils

ENV_TOKEN_FORMAT = "${{{0}}}"
"""
Format string for environment variable tokens used as values
for configuration parameters.
"""
DEFAULT_ENV_DIR = "."
"""
Default directory value used for environment variable tokens
if no value for the environment variable corresponding to the
token is provided.
"""


def get_environment_vars():
    """
    Get map from RiboViz environment variables, to their values,
    if defined. The environment variables are listed in
    :py:const:`riboviz.params.ENV_DIRS`,

    :return: Map from environment variables to values
    :rtype: dict
    """
    return {env: os.environ[env] for env in
            params.ENV_DIRS if env in os.environ}


def update_config_with_env(env_vars, config):
    """
    Replace environment variable tokens with their complementary
    environment variable values in configuration parameter values that
    support environment variable tokens. All other parameters are left
    unchanged.

    :py:const:`riboviz.params.ENV_DIRS` lists parameters whose
    values can include environment variable token
    parameters.

    Tokens are of form ``${<environment_variable>}`` where
    ``<environment_variable>`` is the name of the complementary
    environment variable.

    If any environment variable is undefined then the
    value :py:const.`DEFAULT_ENV_DIR` is inserted.

    :param env_vars: Environment variables and values
    :type env_vars: dict
    :param config: Configuration
    :type config: dict
    """
    env_vars_copy = dict(env_vars)  # Leave original as-is.
    undefined_vars = {env: DEFAULT_ENV_DIR for env in
                      params.ENV_DIRS if env not in env_vars}
    env_vars_copy.update(undefined_vars)
    env_tokens = {ENV_TOKEN_FORMAT.format(env): value
                  for env, value in env_vars_copy.items()}
    for param in params.ENV_PARAMS:
        if param in config and config[param]:
            config[param] = utils.replace_tokens(
                config[param], env_tokens)


def apply_env_to_config(config):
    """
    Replace environment variable tokens with their complementary
    environment variable values in configuration parameter values that
    support environment variable tokens. All other parameters are left
    unchanged.

    See :py:func:`get_environment_vars` and
    :py:func:`update_config_with_env`.

    :param config: Configuration
    :type config: dict
    """
    env_vars = get_environment_vars()
    update_config_with_env(env_vars, config)
