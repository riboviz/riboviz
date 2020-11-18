"""
Python ``logging``-related constants and functions.
"""
import datetime
import logging
import logging.config
import logging.handlers
import os
import yaml

DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), "logging.yaml")
""" Default logging configuration file name. """
LOG_CONFIG_ENV = "RIBOVIZ_LOG_CONFIG"
""" Logging environment variable. """
LOG_FILE = "riboviz.log"
""" Default log file nmae. """


class TimestampedFileHandler(logging.FileHandler):
    """
    Subclass of FileHandler which creates a log file of form
    ``prefixYYYYMMWW-HHMMSSsuffix``.
    """

    def __init__(self, prefix, suffix, mode='a', encoding=None,
                 delay=False):
        """
        Constructor.

        :param prefix: File name prefix
        :type prefix: str or unicode
        :param suffix: File name suffix
        :type suffix: str or unicode
        :param mode: File open mode
        :type mode: str or unicode
        :param encoding: File encoding
        :type encoding: str or unicode
        :param delay: Defer opening file until first output?
        :type delay: bool
        """
        timestamp = str(datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        log_file = prefix + timestamp + suffix
        super(TimestampedFileHandler, self).__init__(
            log_file, mode, encoding, delay)


def configure_logging(config_path=DEFAULT_CONFIG,
                      env_key=LOG_CONFIG_ENV,
                      level=logging.INFO,
                      log_file=LOG_FILE):
    """
    Configure Python logging.

    If an environment variable, whose name is in ``env_key``, is set
    then the configuration file at the location specified by the
    environment variable is used to configure logging.

    If the environment variable is not set, and if ``config_path``, is
    defined, then the configuration file at that location is used to
    configure logging.

    If neither is defined then basic logging is configured, at
    the default level provided, ``level``, into ``log_file``.

    :param config_path: YAML logging configuration file
    :type config_path: str or unicode
    :param env_key: Environment variable name
    :type env_key: str or unicode
    :param level: Default logging level if no configuration file
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
