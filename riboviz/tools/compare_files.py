#!/usr/bin/env python

"""
Command-line utility to compare two files.

Usage: python -m riboviz.tools.compare_files <FILE> <FILE>
"""
import sys
import riboviz.validation

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: compare_files.py <FILE> <FILE>")
        sys.exit(2)
    try:
        riboviz.validation.compare(sys.argv[1], sys.argv[2])
    except AssertionError as e:
        print(e)
        sys.exit(1)
