#!/usr/bin/env python
"""
Upgrade YAML configuration file be compatible with current version of
RiboViz

    usage: upgrade_config_file.py [-h] -i INPUT [-o [OUTPUT]]

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Input YAML configuration file
      -o [OUTPUT], --output [OUTPUT]
                            Output YAML configuration file. If not
                            provided then upgraded content is printed
                            to standard out
"""
import argparse
import sys
from riboviz import upgrade_config_file


def invoke_upgrade_config_file():
    """
    Upgrade YAML configuration file be compatible with current version
    of RiboViz. Parse command-line arguments then invoke
    upgrade_config_file.
    """
    parser = argparse.ArgumentParser(
        description="Upgrade YAML configuration file be compatible with current version of RiboViz")
    parser.add_argument("-i",
                        "--input",
                        dest="input",
                        required=True,
                        help="Input YAML configuration file")
    parser.add_argument("-o",
                        "--output",
                        dest="output",
                        nargs='?',
                        help="Output YAML configuration file. If not provided then upgraded content is printed to standard out")
    options = parser.parse_args()
    input_file = options.input
    output_file = options.output
    upgrade_config_file.upgrade_config_file(input_file, output_file)


if __name__ == "__main__":
    try:
        invoke_upgrade_config_file()
    except Exception as e:
        print(e)
        sys.exit(1)
