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
    file output. This file is used to save having to traverse the
    output SAM file.
* `umi_tools dedup`:
  - Number of reads in output.

The output file is a TSV file with columns:

* `SampleName`
* `Program`
* `File`
* `NumReads`: Number of reads in the file.
* `Description`: Human-readable description of the step.

---

## Example

Here is an example of running `count_reads.py` on workflow files log file produced when running the vignette:

```console
$ python -m riboviz.tools.count_reads \
    -i vignette/output/workflow_files.tsv \
    -o vignette_count_reads.tsv
$ cat vignette_count_reads.tsv
# Created by: RiboViz
# Date: 2020-02-04 01:58:18.166958
# Command-line tool: /home/ubuntu/riboviz/riboviz/tools/count_reads.py
# File: /home/ubuntu/riboviz/riboviz/count_reads.py
# Version: commit 8b3a0f087b982f53f7e88e8dee8ed2b758292b7e date 2020-02-04 01:55:19-08:00
SampleName	Program	File	NumReads	Description
WTnone	input	vignette/input/SRR1042855_s1mi.fastq.gz	963571	Original reads
WTnone	cutadapt	vignette/tmp/WTnone/trim.fq	952343	Reads after removal of sequencing library adapters
WTnone	hisat2	vignette/tmp/WTnone/nonrRNA.fq	466464	rRNA or other contaminating reads removed by alignment to rRNA index files
WTnone	hisat2	vignette/tmp/WTnone/rRNA_map.sam	1430213	Reads with rRNA and other contaminating reads removed by alignment to rRNA index files
WTnone	hisat2	vignette/tmp/WTnone/unaligned.fq	452266	Unaligned reads removed by alignment of remaining reads to ORFs index files
WTnone	hisat2	vignette/tmp/WTnone/orf_map.sam	14516	Reads aligned to ORFs index files
WTnone	riboviz.tools.trim_5p_mismatch	vignette/tmp/WTnone/orf_map_clean.sam	14516	Reads after trimming of 5' mismatches and removal of those with more than 2 mismatches
WT3AT	input	vignette/input/SRR1042864_s1mi.fastq.gz	1374448	Original reads
WT3AT	cutadapt	vignette/tmp/WT3AT/trim.fq	1373362	Reads after removal of sequencing library adapters
WT3AT	hisat2	vignette/tmp/WT3AT/nonrRNA.fq	485226	rRNA or other contaminating reads removed by alignment to rRNA index files
WT3AT	hisat2	vignette/tmp/WT3AT/rRNA_map.sam	2254078	Reads with rRNA and other contaminating reads removed by alignment to rRNA index files
WT3AT	hisat2	vignette/tmp/WT3AT/unaligned.fq	476785	Unaligned reads removed by alignment of remaining reads to ORFs index files
WT3AT	hisat2	vignette/tmp/WT3AT/orf_map.sam	8698e	Reads aligned to ORFs index files
WT3AT	riboviz.tools.trim_5p_mismatch	vignette/tmp/WT3AT/orf_map_clean.sam	8698	Reads after trimming of 5' mismatches and removal of those with more than 2 mismatches
```
