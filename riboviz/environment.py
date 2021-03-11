"""
Environment-related functions.
"""
import os
from riboviz import params


def get_environment_vars():
    """
    Get all RiboViz environment variables, and their values,
    if they are defined.
    The environment variables are:

    * :py:const:`riboviz.params.ENV_RIBOVIZ_SAMPLES`
    * :py:const:`riboviz.params.ENV_RIBOVIZ_ORGANISMS`
    * :py:const:`riboviz.params.ENV_RIBOVIZ_DATA`

    :return: Defined environment variables and values
    :rtype: dict
    """
    return {env: os.environ[env] for env in
            [params.ENV_RIBOVIZ_SAMPLES,
             params.ENV_RIBOVIZ_ORGANISMS,
             params.ENV_RIBOVIZ_DATA] if env in os.environ}
