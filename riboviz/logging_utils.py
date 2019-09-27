"""
Logging-related utilities.
"""
import datetime
import logging
import logging.config
import logging.handlers
import os
import yaml


DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), "logging.yaml")
""" Default logging configuration file. """
LOG_CONFIG_ENV = "RIBOVIZ_LOG_CONFIG"
""" Logging environment variable. """


class TimestampedFileHandler(logging.FileHandler):
    """
    Subclass of FileHandler which creates a log file of form
    prefixYYYYMMWW-HHMMSS-SSSSSSsuffix.
    """
    def __init__(self, prefix, suffix, mode='a', encoding=None, delay=False):
        """
        Constructor.

        :param prefix: file name prefix
        :type prefix: str or unicode
        :param suffix: file name suffix
        :type suffix: str or unicode
        :param mode: mode to open file with
        :type mode: str or unicode
        :param encoding: encoding to use
        :type encoding: str or unicode
        :param delay: defer opening file until first output
        :type delay: bool
        """
        timestamp = str(datetime.datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
        log_file = prefix + timestamp + suffix
        super(TimestampedFileHandler, self).__init__(
            log_file, mode, encoding, delay)


def configure_logging(config_path=DEFAULT_CONFIG,
                      env_key=LOG_CONFIG_ENV,
                      level=logging.INFO,
                      log_file="riboviz.log"):
    """
    Configure logging.

    If an environment variable, env_key, is defined then the
    configuration file at the location specified by the environment
    variable is used to configure logging.

    If it is not set, and if config_path, is defined, then the
    configuration file at that location is used to configure
    logging.

    If neither is defined then basic logging is configured at
    the default level provided, into log_file.

    :param config_path: YAML logging configuration file
    :type config_path: str or unicode
    :param env_key: environment variable with location of YAML logging
    configuration file
    :type env_key: str or unicode
    :param level: default logging level if no configuration file
    :type level: int
    """
    path = config_path
    env_path = os.getenv(env_key, None)
    if env_path is not None:
        path = env_path
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = yaml.load(f, yaml.SafeLoader)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(
            level=level,
            filename=log_file,
            format='%(asctime)s:%(name)s:%(levelname)s: %(message)s')
