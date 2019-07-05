# Developer guide

## Install Python packages for developers

### pylint

Web sites:

* [Pylint](https://www.pylint.org/)
* [BitBucket](https://bitbucket.org/logilab/pylint.org)

Install:
```bash
conda install -y pylint
```

### pycodestyle

Web sites:

* [readthedocs](https://pycodestyle.readthedocs.io/)
* [GitHub](https://github.com/pycqa/pycodestyle)

Install:

```bash
conda install -y pycodestyle
```

### pandas

Web sites:

* [pandas](https://pandas.pydata.org/)
* [GitHub](https://github.com/pandas-dev/pandas)

Install:

```bash
conda install -y pandas
```

### pytest

Web sites:

* [pytest](https://pytest.org/)
* [GitHub](https://github.com/pytest-dev/pytest/)

Install:

```bash
conda install -y pytest
```

---

## Compare files for equality

`pyscripts/compare_files.py` is a script that can compare the files output by any stage of the vignette, `pyscripts/prepRiboViz.py`.

It can be run as follows:

* Python 3:

```bash
python -m pyscripts.compare_files <FILE1> <FILE2>
```

* Python 2 or 3:

```bash
PYTHONPATH=. python pyscripts/compare_files.py <FILE1> <FILE2>
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

---

## Run vignette regression test suite

The vignette (`pyscripts/prepRiboviz.py`) regression test suite compares two directories, each assumed to have `index/` `tmp/` and `output/` directories created by the vignette.

The tests can be run using pytest:

```bash
pytest tests/test_vignette.py --expected=<DIR1> [--actual=<DIR2>]
```

where:

* `--expected=<DIRECTORY>`: directory with expected vignette files.
* `--actual=<DIRECTORY>`: directory to be validated against directory with expected vignette files. Default: `vignette/`

For each file output by the vignette, the files are compared using the same comparisons as for `compare_files.py` above.

Useful pytest flags are:

* `-s`: disable output capture so, for example, `print` messages are shown.
* `-v`: verbose mode, displays names of test functions run.
* `-k`: run a specific test function.

For example:

```bash
pytest -s -k "test_output_bam" tests/test_vignette.py --expected=<DIR1> --actual=<DIR2>
```

---

## Coding style

Regularly run `pylint`, `pycodestyle`, and `2to3` and update the code to adopt the recommendations as far as possible. For example, to run these on `validation.py`:

```bash
pylint riboviz/validation.py
pycodestyle riboviz/validation.py
2to3 riboviz/validation.py
```
