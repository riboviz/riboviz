# count_reads.py reads counter

`riboviz/tools/count_reads.py` is a command-line tool to scan input, temporary and output directories and count the number of reads (sequences) processed by specific stages of a RiboViz workflow. The scan is based on the configuration, directory structure and file patterns used by RiboViz. It outputs a [read counts file](./prep-riboviz-operation.md#read-counts-file).

---

## Usage

```
$ python -m riboviz.tools.count_reads -h

usage: count_reads.py [-h] -c CONFIG_FILE -i INPUT_DIR -t TMP_DIR -o
                      OUTPUT_DIR -r READS_FILE

Scan RiboViz input, temporary and output directories and count the number of
reads (sequences) processed at specific stages of a RiboViz workflow

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        Configuration file
  -i INPUT_DIR, --input-dir INPUT_DIR
                        Input files directory
  -t TMP_DIR, --tmp-dir TMP_DIR
                        Temporary files directory
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        Output files directory
  -r READS_FILE, --reads-file READS_FILE
                        Reads file (output)
```

Inputs:

* `-c CONFIG_FILE`, `--config-file CONFIG_FILE`: Configuration file
* `-i INPUT_DIR`, `--input-dir INPUT_DIR`: Input files directory
* `-t TMP_DIR`, `--tmp-dir TMP_DIR`: Temporary files directory
* `-o OUTPUT_DIR`, `--output OUTPUT_DIR'`: Output files directory

Outputs:

* `-r READS_FILE`, `--reads-file READS_FILE`: Reads file

See [read counts file](./prep-riboviz-operation.md#read-counts-file) for a description of the read counts file format and an example.

---

## Example

Here is an example of running `count_reads.py` on the directories produced when running the vignette:

```console
$ python -m riboviz.tools.count_reads -c vignette/vignette_config.yaml \
    -i vignette/input/ -t vignette/tmp/ -o vignette/output/ \
    -r read_counts.tsv
Created by: RiboViz
Date: 2020-02-12 03:19:04.175143
Command-line tool: /home/ubuntu/riboviz/riboviz/tools/count_reads.py
File: /home/ubuntu/riboviz/riboviz/tools/count_reads.py
Version: commit c1ef2969cb17dbb190f9bbdd9e513539daf38e9f date 2020-02-12 02:11:42-08:00

vignette/input/SRR1042855_s1mi.fastq.gz
vignette/input/SRR1042864_s1mi.fastq.gz
vignette/input/example_missing_file.fastq.gz
[Errno 2] No such file or directory: 'vignette/input/example_missing_file.fastq.gz'
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
```
