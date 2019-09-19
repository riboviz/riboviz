"""
Logging-related utilities.
"""
import datetime
import logging


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
        super(TimestampedFileHandler,self).__init__(
            log_file, mode, encoding, delay)
