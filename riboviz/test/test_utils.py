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


@pytest.mark.parametrize("value", [True, 1, [1], {"a": 1}], ids=str)
@pytest.mark.parametrize("allow_false_empty", [False, True])
def test_value_in_dict(value, allow_false_empty):
    """
    Test that :py:func:`riboviz.utils.value_in_dict` returns `True` if
    a key is present and has one of the given values regardless of the
    value of `allow_false_empty`.

    :param value: Value for key
    :type value: -
    :param allow_false_empty: Value for ``allow_false_empty`` \
    parameter
    :type allow_false_empty: bool
    """
    values = {"A": 1, "B": value, "C": 3}
    assert utils.value_in_dict("B", values, allow_false_empty)


def test_value_in_dict_no_key():
    """
    Test that :py:func:`riboviz.utils.value_in_dict` returns `True` if
    non-existent key always returns ``False``.
    """
    values = {"A": 1, "C": 3}
    assert not utils.value_in_dict("B", values)


@pytest.mark.parametrize("allow_false_empty", [False, True])
def test_value_in_dict_none(allow_false_empty):
    """
    Test that :py:func:`riboviz.utils.value_in_dict` returns ``False`
    if a key has value ``None`` regardless of the value of
    ``allow_false_empty``.

    :param allow_false_empty: Value for ``allow_false_empty`` \
    parameter
    :type allow_false_empty: bool
    """
    values = {"A": 1, "B": None, "C": 3}
    assert not utils.value_in_dict("B", values, allow_false_empty)


@pytest.mark.parametrize("value", [False, [], {}], ids=str)
@pytest.mark.parametrize("allow_false_empty", [False, True])
def test_value_in_dict_allow_false_empty(value, allow_false_empty):
    """
    Test that :py:func:`riboviz.utils.value_in_dict` returns the same
    value as ``allow_false_empty`` if a key is present and has one of
    the given values.

    :param value: Value for key
    :type value: -
    :param allow_false_empty: Value for ``allow_false_empty`` \
    parameter
    :type allow_false_empty: bool
    """
    values = {"A": 1, "B": value, "C": 3}
    is_value_in_dict = utils.value_in_dict("B", values, allow_false_empty)
    assert (not is_value_in_dict or allow_false_empty) and\
           (is_value_in_dict or not allow_false_empty)  # NXOR


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
