# count_reads.py reads counter

`riboviz/tools/count_reads.py` is a command-line tool to process a workflow files log file and count the number of reads (sequences) processed by specific stages of a RiboViz workflow. It outputs a [read counts file](./prep-riboviz-operation.md#read-counts-file).

---

## Usage

```
$ python -m riboviz.tools.count_reads -h
usage: count_reads.py [-h] -i WORKFLOW_FILE -o READS_FILE

Process a workflow files log file and count the number of reads
(sequences) processed at specific stages of a RiboViz workflow

optional arguments:
  -h, --help            show this help message and exit
  -i WORKFLOW_FILE, --input WORKFLOW_FILE
                        Workflow files log file (input)
  -o READS_FILE, --output READS_FILE
                        Reads file (output)
```

Arguments:

* '-h', '--help': show this help message and exit
* '-i WORKFLOW_FILES_LOG_FILE', '--input WORKFLOW_FILES_LOG_FILE':
  workflow files log file (input)
* '-o READS_FILE', '--output READS_FILE': reads file (output)

The input file is assumed to be a TSV file with `riboviz.workflow_files_logger`-consistent columns:

* `SampleName`: Name of the sample to which this file belongs. This is an empty value if the step was not sample-specific (e.g. an input file or a multiplexed FASTQ file).
* `Program`: Program that read/wrote the file. The special token `input` denotes input files.
* `File`: Path to file.
* `Read/Write`: `read` if the file was read, `write` if the file was written.

The files logged in the workflow files log file must exist.

See [read counts file](./prep-riboviz-operation.md#read-counts-file) for a description of the read counts file format and an example.

---

## Example

Here is an example of running `count_reads.py` on workflow files log file produced when running the vignette:

```console
$ python -m riboviz.tools.count_reads \
    -i vignette/output/workflow_files.tsv \
    -o vignette_count_reads.tsv
```
