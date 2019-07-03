# Developer guide

## Compare files

`pyscripts/compare_files.py` is a script that can compare the files output by any stage of the vignette, `pyscripts/prepRiboViz.py`.

It can be run as follows:

* Python 3:

```bash
python -m pyscripts.compare_files <FILE> <FILE>
```

* Python 2 or 3:

```bash
PYTHONPATH=. python pyscripts/compare_files.py <FILE> <FILE>
```

If the files are equivalent an exit code of 0 is returned. If the files are not equivalent an error message is displayed and an exit code of 1 is returned.

For example:

```bash
PYTHONPATH=. python pyscripts/compare_files.py $DIR1/output/WT3AT.h5 $DIR2/output/WT3AT.h5
echo $?
```
```
0
```
```bash
PYTHONPATH=. python pyscripts/compare_files.py $DIR1/output/WT3AT.h5 $DIR2/output/WTnone.h5
```
```
Non-zero return code (1) from h5diff -q .../output/WT3AT.h5 .../output/WTnone.h5
```
```bash
echo $?
```
```
1
```

The following files can be compared:

* `.bam`: compare BAM files for equality, including:
  - Index-specific information:
    - Index statistics.
    - Number of unequal reads without coordinates.
    - Number of mapped alignments.
    - Number of unmapped alignments.
  - Header values for all but "PG".
  - Reference numbers, names and lengths.
  - Reads
  - BAM files are required to have complementary BAI files.
* `.bam.bai`: compare BAI file sizes for equality.
* `.bedgraph`: compare bedGraph file contents for equality.
* `.fq`: compare FASTQ files for equality, including:
  - Both files have the same number of records.
  - All records in file1 are also in file2. The order of records is ignored.
* `.h5`: compare HDF5 files for equality, using `h5diff`.
* `.ht2`: compare HISAT2 file sizes for equality.
* `.pdf`: compare PDF file sizes for equality.
* `.sam`: compare SAM files for equality, including:
  - Header values for all but "PG".
  - Reference numbers, names and lengths.
  - Reads.
* `.tsv`: compare tab-separated (TSV) files for exact equality.
