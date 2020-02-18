"""
Useful functions.
"""
import os
import os.path
import numpy as np
import pandas as pd


def value_in_dict(key, dictionary, allow_false_empty=False):
    """
    Check that a value is in a dictionary and the value is not
    ``None``.

    If dictionary is::

        {
         "A":1,
         "B":None,
         "C":{},"D":[],
         "E":[1],
         "F":True,
         "G":False
        }

    then:

    * ``value_in_dict("A", dictionary)`` is ``True``
    * ``value_in_dict("B", dictionary)`` is ``False``
    * ``value_in_dict("C", dictionary)`` is ``False``
    * ``value_in_dict("D", dictionary)`` is ``False``
    * ``value_in_dict("E", dictionary)`` is ``True``
    * ``value_in_dict("F", dictionary)`` is ``True``
    * ``value_in_dict("G", dictionary)`` is ``False``
    * ``value_in_dict("A", dictionary, True)`` is ``True``
    * ``value_in_dict("B", dictionary, True)`` is ``False``
    * ``value_in_dict("C", dictionary, True)`` is ``True``
    * ``value_in_dict("D", dictionary, True)`` is ``True``
    * ``value_in_dict("E", dictionary, True)`` is ``True``
    * ``value_in_dict("F", dictionary, True)`` is ``True``
    * ``value_in_dict("G", dictionary, True)`` is ``True``

    :param key: Key
    :type key: -
    :param dictionary: Dictionary
    :type dictionary: dict
    :param allow_false_empty: Allow ``False``, empty string, \
    ``list`` or ``dict`` to be considered as an existing value
    :type allow_false_empty: bool
    :return: ``True`` or ``False``
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
    first ``.`` in the file name. For example, given
    ``example.fastq.gz`` return ``fastq.gz``, given  ``example.fastq``
    return ``fastq``, given ``example`` return ``''``. The extension
    is returned in lower-case.

    :param file_name: File name
    :type file_name: str or uniecode
    :return: Extension
    :rtype: str or unicode
    """
    file_type = ".".join(os.path.basename(file_name).split(".")[1:])
    return file_type.lower()


def equal_file_names(file1, file2):
    """
    Compare local names of two files each of which must exist and be a
    file.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: If file do not exist, are not files or
    their names differ
    """
    local_file1 = os.path.split(file1)[1].lower()
    local_file2 = os.path.split(file2)[1].lower()
    assert os.path.exists(file1) and os.path.isfile(file1),\
        "File %s does not exist or is not a file"
    assert os.path.exists(file2) and os.path.isfile(file2),\
        "File %s does not exist or is not a file"
    assert local_file1 == local_file2,\
        "Unequal file names: %s, %s" % (local_file1, local_file2)


def equal_file_sizes(file1, file2):
    """
    Compare sizes of two files.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :raise AssertionError: If the file sizes differ
    :raise Exception: If problems arise when accessing the files
    """
    stat1 = os.stat(file1)
    stat2 = os.stat(file2)
    assert stat1.st_size == stat2.st_size,\
        "Unequal file sizes: %s, %s" % (file1, file2)


def equal_dataframes(data1, data2, tolerance=0.0001):
    """
    Compare two Pandas data frames for equality. The data frames are
    expected to be two dimensional i.e. rows and columns.

    The data frames are compared column-by-column:

    * ``float64`` columns are converted to numpy arrays then tested for
      equality to within the given tolerance using
      ``numpy.allclose``. This is used instead of
      ``pandas.testing.assert_frame_equal`` as there is an issue with
      how that function handles precision (see
      'pandas.testing.assert_frame_equal doesn't do precision
      according to the doc' #25068,
      https://github.com/pandas-dev/pandas/issues/25068). In addition,
      ``NAN`` values are considered to be equal.
    * All other columns (``object``, ``int64``, ``bool``,
      ``datetime64``, ``timedelta``) are compared for exact equality
      using ``pandas.core.series.Series.equals``.

    :param data1: dataframe
    :type data1: pandas.core.frame.DataFrame
    :param data2: dataframe
    :type data2: pandas.core.frame.DataFrame
    :param tolerance: Tolerance for floating point comparisons
    :type tolerance: float
    :raise AssertionError: If the data frames differ in their content
    """
    assert data1.shape == data2.shape,\
        "Unequal shape: %s, %s"\
        % (str(data1.shape), str(data2.shape))
    assert data1.columns.equals(data2.columns),\
        "Unequal column names: %s, %s"\
        % (str(data1.columns), str(data2.columns))
    for column in data1.columns:
        column1 = data1[column]
        column2 = data2[column]
        if column1.dtype in (int, float) and column2.dtype in (int, float):
            column_data1 = column1.to_numpy()
            column_data2 = column2.to_numpy()
            assert np.allclose(column_data1,
                               column_data2,
                               rtol=0,
                               atol=tolerance,
                               equal_nan=True),\
                "Unequal column values: %s" % column
        else:
            assert column1.equals(column2),\
                "Unequal column values: %s" % column


def equal_tsv(file1, file2, tolerance=0.0001, comment="#"):
    """
    Compare two tab-separated (TSV) files for equality. This function
    uses :py:func:`equal_dataframes`.

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :param tolerance: Tolerance for floating point comparisons
    :type tolerance: float
    :param comment: Comment prefix
    :type comment: str or unicode
    :raise AssertionError: If files differ in their contents
    :raise Exception: If problems arise when loading the files
    """
    data1 = pd.read_csv(file1, sep="\t", comment=comment)
    data2 = pd.read_csv(file2, sep="\t", comment=comment)
    try:
        equal_dataframes(data1, data2, tolerance)
    except AssertionError as error:
        # Add file names to error message.
        message = error.args[0]
        message += " in file: " + str(file1) + ":" + str(file2)
        error.args = (message,)
        raise
