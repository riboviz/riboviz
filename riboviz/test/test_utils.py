"""
utils tests.
"""
import pytest
from riboviz.utils import list_to_str
from riboviz.utils import value_in_dict


def test_list_to_str_empty():
    """
    Test list_to_str with an empty list.
    """
    list_str = list_to_str([])
    assert list_str == "", "Unexpected string {}".format(list_str)


def test_list_to_str():
    """
    Test list_to_str with a list of str.
    """
    list_str = list_to_str(["ab", "cd", "ef"])
    assert list_str == "ab cd ef", "Unexpected string {}".format(list_str)


def test_list_to_str_int():
    """
    Test list_to_str with a list of int.
    """
    list_str = list_to_str([1, 2, 3])
    assert list_str == "1 2 3", "Unexpected string {}".format(list_str)


@pytest.mark.parametrize("value", [True, 1, [1], {"a": 1}])
@pytest.mark.parametrize("allow_false_empty", [False, True])
def test_value_in_dict(value, allow_false_empty):
    """
    Test value_in_dict returns True if a key is present and has one of
    the given values regardless of the value of allow_false_empty.

    :param value: value for key
    :type value: -
    :param allow_false_empty: Value for allow_false_empty parameter of
    value_in_dict
    :type allow_false_empty: bool
    """
    values = {"A": 1, "B": value, "C": 3}
    assert value_in_dict("B", values, allow_false_empty)


def test_value_in_dict_no_key():
    """
    Test value_in_dict with a non-existent key always returns False.
    """
    values = {"A": 1, "C": 3}
    assert not value_in_dict("B", values)


@pytest.mark.parametrize("allow_false_empty", [False, True])
def test_value_in_dict_none(allow_false_empty):
    """
    Test value_in_dict returns False if a key has value None
    regardless of the value of allow_false_empty.

    :param allow_false_empty: Value for allow_false_empty parameter of
    value_in_dict
    :type allow_false_empty: bool
    """
    values = {"A": 1, "B": None, "C": 3}
    assert not value_in_dict("B", values, allow_false_empty)


@pytest.mark.parametrize("value", [False, [], {}])
@pytest.mark.parametrize("allow_false_empty", [False, True])
def test_value_in_dict_allow_false_empty(value, allow_false_empty):
    """
    Test value_in_dict returns the same value as allow_false_empty if
    a key is present and has one of the given values.

    :param value: value for key
    :type value: -
    :param allow_false_empty: Value for allow_false_empty parameter of
    value_in_dict
    :type allow_false_empty: bool
    """
    values = {"A": 1, "B": value, "C": 3}
    is_value_in_dict = value_in_dict("B", values, allow_false_empty)
    assert (not is_value_in_dict or allow_false_empty) and\
           (is_value_in_dict or not allow_false_empty)  # NXOR
