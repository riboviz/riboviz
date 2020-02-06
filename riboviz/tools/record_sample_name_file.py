"""
Record the association between a sample name and a sample file.

Usage:

    record_sample_name_file.py [-h] -r RECORD_FILE -n NAME -f FILE

* '-r RECORD_FILE', '--record-file RECORD_FILE': Record file
* '-n NAME', '--name NAME':  Sample name
* '-f FILE', '--file FILE':  Sample file

The record file is a TSV file with columns:

* 'SampleName': Name of the sample to which this file belongs.
* 'File': Path to file.

If the file already exists the record will be appended to the file.
"""
import argparse
from riboviz import sample_names_files
from riboviz import provenance


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Record the association between a sample name and a sample file.")
    parser.add_argument("-r",
                        "--record-file",
                        dest="record_file",
                        required=True,
                        help="Record file")
    parser.add_argument("-n",
                        "--name",
                        dest="sample_name",
                        required=True,
                        help="Sample name")
    parser.add_argument("-f",
                        "--file",
                        dest="sample_file",
                        required=True,
                        help="Sample file")
    options = parser.parse_args()
    return options


def invoke_record_sample():
    """
    Parse command-line options then invoke "record_sample".
    """
    print(provenance.get_provenance_str(__file__))
    options = parse_command_line_options()
    record_file = options.record_file
    sample_name = options.sample_name
    sample_file = options.sample_file
    sample_names_files.record_sample(record_file, sample_name, sample_file)


if __name__ == "__main__":
    invoke_record_sample()
