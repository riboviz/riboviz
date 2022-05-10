"""
:py:mod:`riboviz.utils` tests.
"""
import pytest
from riboviz import utils


def test_list_to_str_empty():
    """
    Test :py:func:`riboviz.utils.list_to_str` with an empty list.
    """
    list_str = utils.list_to_str([])
    assert list_str == "", "Unexpected string {}".format(list_str)


def test_list_to_str():
    """
    Test :py:func:`riboviz.utils.list_to_str` with a list of ``str``.
    """
    list_str = utils.list_to_str(["ab", "cd", "ef"])
    assert list_str == "ab cd ef", "Unexpected string {}".format(list_str)


def test_list_to_str_int():
    """
    Test :py:func:`riboviz.utils.list_to_str` with a list of ``int``.
    """
    list_str = utils.list_to_str([1, 2, 3])
    assert list_str == "1 2 3", "Unexpected string {}".format(list_str)


def get_file_ext_name_dot_ext_ext():
    """
    Test :py:func:`riboviz.utils.get_file_ext` with
    ``example.fastq.gz``.
    """
    assert utils.get_file_ext("example.fastq.gz") == "fastq.gz"


def get_file_ext_name_dot_ext():
    """
    Test :py:func:`riboviz.utils.get_file_ext` with
    ``example.fastq``.
    """
    assert utils.get_file_ext("example.fastq") == "fastq"


def get_file_ext_name():
    """
    Test :py:func:`riboviz.utils.get_file_ext` with
    ``example``.
     """
    assert utils.get_file_ext("example") == ""


@pytest.mark.parametrize(
    "string,tokens,expected",
    [("", {}, ""),
     ("abcdef", {}, "abcdef"),
     ("", {"ab": "GH", "ef": "IJ"}, ""),
     ("abcdef", {"ab": "GH", "ef": "IJ"}, "GHcdIJ"),
     ("abcdef", {"pq": "QP", "ef": "IJ"}, "abcdIJ"),
     ("abcdef", {"pq": "QP", "rs": "SR"}, "abcdef")],
    ids=["empty-string-empty-tokens",
         "string-empty-tokens",
         "empty-string-tokens",
         "string-all-token",
         "string-some-tokens",
         "string-no-tokens"])
def test_replace_tokens(string, tokens, expected):
    """
    Test :py:func:`riboviz.utils.replace_tokens`.

    :param string: String
    :type string: str or unicode
    :param tokens: Map from tokens to substrings
    :type tokens: dict
    :param expected: Expected result
    :type expected: str or unicode
     """
    assert utils.replace_tokens(string, tokens) == expected
