"""
Utilities.
"""
import os


def value_in_dict(key, dictionary, allow_false_empty=False):
    """
    Check that a value is in a dictionary and the value is not None.

    If dictionary is

        {"A":1,"B":None,"C":{},"D":[],"E":[1],"F":True,"G":False}

    then

        value_in_dict("A", dictionary): True
        value_in_dict("B", dictionary): False
        value_in_dict("C", dictionary): False
        value_in_dict("D", dictionary): False
        value_in_dict("E", dictionary): True
        value_in_dict("F", dictionary): True
        value_in_dict("G", dictionary): False

        value_in_dict("A", dictionary, True): True
        value_in_dict("B", dictionary, True): False
        value_in_dict("C", dictionary, True): True
        value_in_dict("D", dictionary, True): True
        value_in_dict("E", dictionary, True): True
        value_in_dict("F", dictionary, True): True
        value_in_dict("G", dictionary, True): True

    :param key: Key
    :type key: -
    :param dictionary: Dictionary
    :type dictionary: dict
    :param allow_false_empty: Consider False, empty string, list, dict
    to be existant
    :type allow_false_empty: bool
    :return: True or False
    :rtype: bool
    """
    is_in = key in dictionary and dictionary[key] is not None
    if not allow_false_empty:
        is_in = is_in and bool(dictionary[key])
    return is_in


def list_to_str(lst):
    """
    Convert list to space-delimited string.

    :param lst: list
    :type lst: list
    :return: list as string
    :rtype: str or unicode
    """
    return ' '.join(map(str, lst))


def get_file_ext(file_name):
    """
    Given a file name return full file extension, everything after the
    first "." in the file name. For example, for 'example.fastq.gz'
    return 'fastq.gz', for 'example.fastq' return 'fastq', for 'example'
    return ''. The extension is returned in lower-case.

    :param file_name: File name
    :type file_name: str or unicode
    :return: extension
    :rtype: str or unicode
    """
    file_type = ".".join(os.path.basename(file_name).split(".")[1:])
    return file_type.lower()
