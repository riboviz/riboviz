"""
Environment-related functions.
"""
import os
from riboviz import params


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


def get_environment_tokens():
    """
    Get all RiboViz environment variable tokens, to the
    values of the corresponding environment variables,
    if defined. The environment variables are:

    * :py:const:`riboviz.params.ENV_RIBOVIZ_SAMPLES`
    * :py:const:`riboviz.params.ENV_RIBOVIZ_ORGANISMS`
    * :py:const:`riboviz.params.ENV_RIBOVIZ_DATA`

    The tokens are of form ``${ENVIRONMENT_VARIABLE}``.

    :return: Map from environment variable tokens to values
    :rtype: dict
    """
    env_values = get_environment_vars()
    return {"${{{0}}}".format(env): value
            for env, value in env_values.items()}
