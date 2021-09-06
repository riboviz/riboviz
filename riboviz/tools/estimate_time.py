#!/usr/bin/env python
"""
Estimate the execution time of one tasks, based on the execution
history of some down sampled data.

Usage::

    python -m riboviz.tools.estimate_time [-h]
        -s SAMPLE_SIZE -g GFF_LINES -p NUM_PROCESS -t TASK

    -h, --help            show this help message and exit
    -s SAMPLE_SIZE, --sample-size SAMPLE_SIZE
                          Sample File (*.fastq.gz) size in bytes
    -g GFF_LINES, --gff-lines GFF_LINES
                          Lines of the gff (*.gff3) file (including comments)
    -p NUM_PROCESS, --num-process NUM_PROCESS
                          The number of process that the task will run with.
                          Must be an integer between 1 and 32.
    -t TASK, --task TASK
                          The task name, currently only bamtoh5, cutadapt,
                          hisat2orf, hisat2rrna and samviewsort are supported.
                          Case insensitive.
    Return value: The estimate time in seconds. If the input is illegal,
                  2 hour (7200) will be returned.

Example::

    $ python -m riboviz.tools.estimate_time \
      -s 1000000000 -g 15000 -p 8 -t bamToH5

"""
import argparse
from riboviz import estimate_time, provenance


def parse_command_line_options():
    """
    Parse command-line options.

    :returns: command-line options
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Estimate the execution time of some tasks in the workflow")
    parser.add_argument("-s",
                        "--sample-size",
                        dest="sample_size",
                        required=False,
                        help="Size of the sample file")
    parser.add_argument("-g",
                        "--gff-lines",
                        dest="gff_lines",
                        required=False,
                        help="The lines of gff file")
    parser.add_argument("-p",
                        "--num_process",
                        dest="num_process",
                        required=False,
                        help="Number of process that the task will use")
    parser.add_argument("-t",
                        "--task",
                        dest="task",
                        required=False,
                        help="Task name, case insensitive")
    options = parser.parse_args()
    return options


def invoke_count_reads():
    """
    Parse command-line options then invoke
    :py:func:`riboviz.count_reads.count_reads`.
    """
    print(provenance.write_provenance_to_str(__file__))
    options = parse_command_line_options()
    try:
        sample_size = options.sample_size
        gff_lines = options.gff_lines
        num_process = options.num_process
        task = options.task.lower()
        time = estimate_time.estimate_time(float(sample_size), int(gff_lines), int(num_process), task) / 60
        if time > 126 or time < 0:
            time = 126
    except:
        exit(120)
    exit(int(time+1))


if __name__ == "__main__":
    invoke_count_reads()
