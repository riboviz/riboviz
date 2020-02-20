# `riboviz.tools.compare_files` file comparison tool

`riboviz.tools.compare_files` is a command-line tool that can compare the files output by any stage of the workflow, `riboviz.tools.prep_riboviz`.

---

## Usage:

```console
$ python -m riboviz.tools.compare_files -i INPUT_FILE -o OUTPUT_FILE [-n]
```

where:

* `-1 FILE1`, `--file1 FILE1`: File1
* `-2 FILE2`, `--file2 FILE2`: File2
* `-n`, `--names`: Compare file names

If the files are equivalent an exit code of 0 is returned. If the files are not equivalent an error message is displayed and an exit code of 1 is returned.

For example:

```console
$ python -m riboviz.tools.compare_files -i $DIR1/output/WT3AT.h5 -o $DIR2/output/WT3AT.h5
$ echo $?
0

$ python -m riboviz.tools.compare_files -i $DIR1/output/WT3AT.h5 -o $DIR2/output/WTnone.h5
AssertionError: Non-zero return code (1) from h5diff -q $DIR1/output/WT3AT/WT3AT.h5 $DIR2/output/WTnone/WTnone.h5
$ echo $?
1
```

The following files can be compared:

* `.bam`: compare sorted BAM files for equality, including:
  - Index-specific information:
    - Index statistics.
    - Number of unequal reads without coordinates.
    - Number of mapped alignments.
    - Number of unmapped alignments.
  - Category, version, compression, description.
  - Header values for all but "PG".
  - Reference numbers, names and lengths.
  - Reads
  - BAM files are required to have complementary BAI files.
  - BAM files are expected to be sorted by leftmost coordinate position.
* `.bam.bai`: compare BAI file sizes for equality.
* `.bedgraph`: compare bedGraph file contents for equality.
* `.fq`: compare FASTQ files for equality, including:
  - Both files have the same number of records.
  - All records in file1 are also in file2. The order of records is ignored.
* `.h5`: compare HDF5 files for equality, using `h5diff`.
* `.ht2`: compare HISAT2 file sizes for equality.
* `.pdf`: compare PDF file sizes for equality.
* `.sam`: compare sorted SAM files for equality, including:
  - Category, version, compression, description.
  - Header values for all but "PG".
  - Reference numbers, names and lengths.
  - Reads.
  - SAM files are expected to be sorted by leftmost coordinate position.
* `.tsv`: compare tab-separated (TSV) files for exact equality.
