# count_reads.py reads counter

`riboviz/tools/count_reads.py` is a command-line tool to process a workflow files log file and count the number of reads (sequences) processed by specific stages of a RiboViz workflow.

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

The input file is assumed to be a TSV file with
`riboviz.workflow_files_logger`-consistent columns:

* `SampleName`: Name of the sample to which this file belongs. This is
  an empty value if the step was not sample-specific (e.g. creating
  index files or demultiplexing a multiplexed FASTQ file).
* `Program`: Program that read/wrote the file. The special token
  `input` denotes input files.
* `File`: Path to file read/written.
* `Read/Write`: `read` if the file was read, `write` if the file was
  written.

The files logged in the workflow files log file must exist.

The following information is extracted:

* Input files:
  - Number of reads in sample FASTQ files used as inputs (if
    non-multiplexed samples are used).
  - Number of reads in multiplexed FASTQ file (if multiplexed
    samples are used).
* `cutadapt`:
  - Number of reads in FASTQ file output as recorded in the FASTQ
    file output.
* `riboviz.tools.demultiplex_fastq`
  - Number of reads in demultiplexed FASTQ files and unassigned reads,
    as recorded in the `num_reads.tsv` file output. This file is used
    to save having to traverse each of the output FASTQ files.
* `hisat2`:
  - Number of reads in SAM file and FASTQ files output.
* `riboviz.tools.trim_5p_mismatch`:
  - Number of reads in SAM file output as recorded in the TSV summary
    file output.This file is used to save having to traverse each of
    the output SAM files.
* `umi_tools dedup`:
  - Number of reads in output.

The output file is a TSV file with columns:

* `SampleName`
* `Program`
* `File`
* `NumReads`: Number of reads in the file.
* `Description`: Human-readable description of the step.
