#!/usr/bin/env python
"""
Scan input, temporary and output directories and count the number of
reads (sequences) processed by specific stages of a workflow.

Usage::

    python -m riboviz.tools.count_reads [-h]
        -c CONFIG_FILE -i INPUT_DIR -t TMP_DIR -o OUTPUT_DIR
        -r READS_FILE

    -h, --help            show this help message and exit
    -c CONFIG_FILE, --config-file CONFIG_FILE
                          Configuration file
    -i INPUT_DIR, --input-dir INPUT_DIR
                          Input directory
    -t TMP_DIR, --tmp-dir TMP_DIR
                          Temporary directory
    -o OUTPUT_DIR, --output OUTPUT_DIR
                          Output directory
    -r READS_FILE, --reads-file READS_FILE
                          Reads file (output)

Example::

    $ python -m riboviz.tools.count_reads
      -c vignette/vignette_config.yaml
      -i vignette/input/ -t vignette/tmp/ -o vignette/output/
      -r read_counts.tsv
    ...
    vignette/input/SRR1042855_s1mi.fastq.gz
    vignette/input/SRR1042864_s1mi.fastq.gz
    vignette/input/example_missing_file.fastq.gz
    [Errno 2] No such file or directory:
        'vignette/input/example_missing_file.fastq.gz'
    vignette/tmp/WT3AT/trim.fq
    vignette/tmp/WT3AT/nonrRNA.fq
    vignette/tmp/WT3AT/rRNA_map.sam
    vignette/tmp/WT3AT/unaligned.fq
    vignette/tmp/WT3AT/orf_map.sam
    vignette/tmp/WT3AT/trim_5p_mismatch.tsv
    vignette/tmp/WTnone/trim.fq
    vignette/tmp/WTnone/nonrRNA.fq
    vignette/tmp/WTnone/rRNA_map.sam
    vignette/tmp/WTnone/unaligned.fq
    vignette/tmp/WTnone/orf_map.sam
    vignette/tmp/WTnone/trim_5p_mismatch.tsv

See :py:func:`riboviz.count_reads.count_reads` for information on what
files are read and how the reads are counted.
"""
import argparse
from riboviz import count_reads
from riboviz import provenance


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Scan input, temporary and output directories and count the number of reads (sequences) processed by specific stages of a workflow")
    parser.add_argument("-c",
                        "--config-file",
                        dest="config_file",
                        required=True,
                        help="Configuration file")
    parser.add_argument("-i",
                        "--input-dir",
                        dest="input_dir",
                        required=True,
                        help="Input directory")
    parser.add_argument("-t",
                        "--tmp-dir",
                        dest="tmp_dir",
                        required=True,
                        help="Temporary directory")
    parser.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="Output directory")
    parser.add_argument("-r",
                        "--reads-file",
                        dest="reads_file",
                        required=True,
                        help="Reads file (output)")
    options = parser.parse_args()
    return options


def main():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.count_reads.count_reads`.
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    config_file = options.config_file
    input_dir = options.input_dir
    tmp_dir = options.tmp_dir
    output_dir = options.output_dir
    reads_file = options.reads_file
    count_reads.count_reads(
        config_file, input_dir, tmp_dir, output_dir, reads_file)


if __name__ == "__main__":
    main()
