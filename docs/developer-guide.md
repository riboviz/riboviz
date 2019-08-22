# Developer guide

## Install Python packages for developers

### pylint

Web sites:

* [Pylint](https://www.pylint.org/)
* [BitBucket](https://bitbucket.org/logilab/pylint.org)

Install:
```console
$ conda install -y pylint
```

### pycodestyle

Web sites:

* [readthedocs](https://pycodestyle.readthedocs.io/)
* [GitHub](https://github.com/pycqa/pycodestyle)

Install:

```console
$ conda install -y pycodestyle
```

### pandas

Web sites:

* [pandas](https://pandas.pydata.org/)
* [GitHub](https://github.com/pandas-dev/pandas)

Install:

```console
$ conda install -y pandas
```

### pytest

Web sites:

* [pytest](https://pytest.org/)
* [GitHub](https://github.com/pytest-dev/pytest/)

Install:

```console
$ conda install -y pytest
```

---

## Branching model

The `master` branch is a stable branch. Releases are created via tags from the master branch.

The `develop` branch is for new functionality. At regular intervals this will be merged into the `master` branch.

For branches relating to issues (i.e. new features, significant refactorings or bug fixes), the naming scheme `<name>-<#>` is used, where:

* `<name>` is a short name for the issue. This should be lower-case, with `-` used as a delimiter if desired.
* `<#>` is the number of the issue as in the GitHub issue tracker.

For example, for the issue "Investigate cutadapt -j flag #43" the branch name was `configure-cutadapt-cores-43`.

To request that a branch be merged into the `develop` branch create a new pull request.

---

## Compare files for equality

`pyscripts/compare_files.py` is a script that can compare the files output by any stage of the vignette, `pyscripts/prep_riboviz.py`.

It can be run as follows:

* Python 3:

```console
$ python -m pyscripts.compare_files <FILE1> <FILE2>
```

* Python 2 or 3:

```console
$ PYTHONPATH=. python pyscripts/compare_files.py <FILE1> <FILE2>
```

If the files are equivalent an exit code of 0 is returned. If the files are not equivalent an error message is displayed and an exit code of 1 is returned.

For example:

```console
$ PYTHONPATH=. python pyscripts/compare_files.py $DIR1/output/WT3AT.h5 $DIR2/output/WT3AT.h5
$ echo $?
0

$ PYTHONPATH=. python pyscripts/compare_files.py $DIR1/output/WT3AT.h5 $DIR2/output/WTnone.h5
Non-zero return code (1) from h5diff -q .../output/WT3AT.h5 .../output/WTnone.h5
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

---

## Run vignette regression test suite

The vignette (`pyscripts/prep_riboviz.py`) regression test suite compares two directories, each assumed to have `index/` `tmp/` and `output/` directories created by the vignette.

The tests can be run using pytest:

```console
$ pytest tests/test_vignette.py --expected=<DIR1> [--actual=<DIR2>]
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

```console
$ pytest -s -k "test_output_bam" tests/test_vignette.py --expected=<DIR1> --actual=<DIR2>
```

A complete run looks like the following

```console
$ pytest -v tests/test_vignette.py --expected=$HOME/expected-vignette --actual=$HOME/actual-vignette
============================= test session starts ==============================
...
collected 51 items                                                             

tests/test_vignette.py::test_index[YAL_CDS_w_250-1] PASSED               [  1%]
tests/test_vignette.py::test_index[YAL_CDS_w_250-2] PASSED               [  3%]
tests/test_vignette.py::test_index[YAL_CDS_w_250-3] PASSED               [  5%]
tests/test_vignette.py::test_index[YAL_CDS_w_250-4] PASSED               [  7%]
tests/test_vignette.py::test_index[YAL_CDS_w_250-5] PASSED               [  9%]
tests/test_vignette.py::test_index[YAL_CDS_w_250-6] PASSED               [ 11%]
tests/test_vignette.py::test_index[YAL_CDS_w_250-7] PASSED               [ 13%]
tests/test_vignette.py::test_index[YAL_CDS_w_250-8] PASSED               [ 15%]
tests/test_vignette.py::test_index[yeast_rRNA-1] PASSED                  [ 17%]
tests/test_vignette.py::test_index[yeast_rRNA-2] PASSED                  [ 19%]
tests/test_vignette.py::test_index[yeast_rRNA-3] PASSED                  [ 21%]
tests/test_vignette.py::test_index[yeast_rRNA-4] PASSED                  [ 23%]
tests/test_vignette.py::test_index[yeast_rRNA-5] PASSED                  [ 25%]
tests/test_vignette.py::test_index[yeast_rRNA-6] PASSED                  [ 27%]
tests/test_vignette.py::test_index[yeast_rRNA-7] PASSED                  [ 29%]
tests/test_vignette.py::test_index[yeast_rRNA-8] PASSED                  [ 31%]
tests/test_vignette.py::test_tmp_fq[WT3AT-nonrRNA] PASSED                [ 33%]
tests/test_vignette.py::test_tmp_fq[WT3AT-trim] PASSED                   [ 35%]
tests/test_vignette.py::test_tmp_fq[WT3AT-unaligned] PASSED              [ 37%]
tests/test_vignette.py::test_tmp_fq[WTnone-nonrRNA] PASSED               [ 39%]
tests/test_vignette.py::test_tmp_fq[WTnone-trim] PASSED                  [ 41%]
tests/test_vignette.py::test_tmp_fq[WTnone-unaligned] PASSED             [ 43%]
tests/test_vignette.py::test_tmp_sam[WT3AT-orf_map_clean] PASSED         [ 45%]
tests/test_vignette.py::test_tmp_sam[WT3AT-orf_map] PASSED               [ 47%]
tests/test_vignette.py::test_tmp_sam[WT3AT-rRNA_map] PASSED              [ 49%]
tests/test_vignette.py::test_tmp_sam[WTnone-orf_map_clean] PASSED        [ 50%]
tests/test_vignette.py::test_tmp_sam[WTnone-orf_map] PASSED              [ 52%]
tests/test_vignette.py::test_tmp_sam[WTnone-rRNA_map] PASSED             [ 54%]
tests/test_vignette.py::test_output_bai[WT3AT] PASSED                    [ 56%]
tests/test_vignette.py::test_output_bai[WTnone] PASSED                   [ 58%]
tests/test_vignette.py::test_output_bam[WT3AT] PASSED                    [ 60%]
tests/test_vignette.py::test_output_bam[WTnone] PASSED                   [ 62%]
tests/test_vignette.py::test_output_bedgraph[WT3AT-minus] PASSED         [ 64%]
tests/test_vignette.py::test_output_bedgraph[WT3AT-plus] PASSED          [ 66%]
tests/test_vignette.py::test_output_bedgraph[WTnone-minus] PASSED        [ 68%]
tests/test_vignette.py::test_output_bedgraph[WTnone-plus] PASSED         [ 70%]
tests/test_vignette.py::test_output_h5[WT3AT] PASSED                     [ 72%]
tests/test_vignette.py::test_output_h5[WTnone] PASSED                    [ 74%]
tests/test_vignette.py::test_output_tsv[WT3AT-3nt_periodicity] PASSED    [ 76%]
tests/test_vignette.py::test_output_tsv[WT3AT-codon_ribodens] PASSED     [ 78%]
tests/test_vignette.py::test_output_tsv[WT3AT-pos_sp_nt_freq] PASSED     [ 80%]
tests/test_vignette.py::test_output_tsv[WT3AT-pos_sp_rpf_norm_reads] PASSED [ 82%]
tests/test_vignette.py::test_output_tsv[WT3AT-read_lengths] PASSED       [ 84%]
tests/test_vignette.py::test_output_tsv[WT3AT-tpms] PASSED               [ 86%]
tests/test_vignette.py::test_output_tsv[WTnone-3nt_periodicity] PASSED   [ 88%]
tests/test_vignette.py::test_output_tsv[WTnone-codon_ribodens] PASSED    [ 90%]
tests/test_vignette.py::test_output_tsv[WTnone-pos_sp_nt_freq] PASSED    [ 92%]
tests/test_vignette.py::test_output_tsv[WTnone-pos_sp_rpf_norm_reads] PASSED [ 94%]
tests/test_vignette.py::test_output_tsv[WTnone-read_lengths] PASSED      [ 96%]
tests/test_vignette.py::test_output_tsv[WTnone-tpms] PASSED              [ 98%]
tests/test_vignette.py::test_output_tpms_collated_tsv PASSED             [100%]

========================= 51 passed in 614.86 seconds ==========================
```

---

## Coding style

Regularly run `pylint`, `pycodestyle`, and `2to3` and update the code to adopt the recommendations as far as possible. For example, to run these on `validation.py`:

```console
$ pylint riboviz/validation.py
$ pycodestyle riboviz/validation.py
$ 2to3 riboviz/validation.py
```

---

## Repository structure

```
data/      # Data files used by scripts and vignette
pyscripts/ # Python scripts invoked by vignette and vignette itself
riboviz/   # Python package code
rscripts/  # R scripts invoked by vignette
rmarkdown/ # Rmarkdown scripts for data preprocessing
tests/     # Python tests including vignette regression test
website/   # RiboViz Shiny server code and data
```

---

## Data file origins

The following data files were created either manually or via the use of scripts outwith or within the repository.

```
data/yeast_CDS_w_250utrs.fa
data/yeast_CDS_w_250utrs.gff3
data/yeast_codon_pos_i200.RData
```

Created by a run of [script_for_transcript_annotation.Rmd](../rmarkdown/script_for_transcript_annotation.Rmd) on third-party data. See [Inputs](./run-vignette.md#inputs) in [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md) for details. These files are used as inputs to RiboViz.

```
vignette/input/yeast_YAL_CDS_w_250utrs.fa
vignette/input/yeast_YAL_CDS_w_250utrs.gff3
vignette/input/yeast_rRNA_R64-1-1.fa
vignette/input/SRR1042855_s1mi.fastq.gz
vignette/input/SRR1042864_s1mi.fastq.gz
```

Created manually from third-party data or the foregoing `data/` FILES. See [Inputs](./run-vignette.md#inputs) in [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md) for details. These files are used as inputs to RiboViz.

```
data/testdata_trim_5p_mismatch.sam
data/testdata_trim_5pos5neg.sam
```

Created by running RiboViz on the data in `vignette/`, copying and pasting lines from SAM files produced, then manually editing the lines to produce the desired range of outcomes. These files are used for testing [trim_5p_mismatch.py](../pyscripts/trim_5p_mismatch.py).
