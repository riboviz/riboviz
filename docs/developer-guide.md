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

## Install R packages for developers

### lintr

Web sites:
* [lintr package on CRAN](https://cran.r-project.org/package=lintr)
* [GitHub](https://github.com/jimhester/lintr)

Install:

```console
$ R
> install.packages("lintr")
# load the package before use with:
 # library("lintr")
```

### styleR

Web sites:
* [StyleR package documentation](https://styler.r-lib.org/)
* [GitHub](https://github.com/r-lib/styler)

Install:

```console
$ R
> install.packages("styler")
# load the package before use with:
 # library("styler")
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

`riboviz/tools/compare_files.py` is a script that can compare the files output by any stage of the workflow, `riboviz/tools/prep_riboviz.py`.

It can be run as follows:

* Either:

```console
$ python -m riboviz.tools.compare_files <FILE1> <FILE2>
```

* Or:

```console
$ PYTHONPATH=. python riboviz/tools/compare_files.py <FILE1> <FILE2>
```

If the files are equivalent an exit code of 0 is returned. If the files are not equivalent an error message is displayed and an exit code of 1 is returned.

For example:

```console
$ PYTHONPATH=. python riboviz/tools/compare_files.py $DIR1/output/WT3AT.h5 $DIR2/output/WT3AT.h5
$ echo $?
0

$ PYTHONPATH=. python riboviz/tools/compare_files.py $DIR1/output/WT3AT.h5 $DIR2/output/WTnone.h5
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

The vignette regression test suite optionally runs `prep_riboviz.py` on the vignette data, in `vignette/`, then compares the results, in `vignette/`, to a directory of pre-calculated results.

The tests can be run using pytest:

```console
$ pytest riboviz/test/regression/test_vignette.py \
        --expected=<DIRECTORY> \
        [--skip-workflow]
```

where:

* `--expected=<DIRECTORY>`: directory with expected vignette files. This is assumed to have `index/` `tmp/` and `output/` directories.
* `--skip-workflow`: request that the `prep_riboviz.py` workflow **not** be run, instead use existing data files in `vignette/` for testing.

For each file output by the vignette, the files are compared using the same comparisons as for `compare_files.py` above.

Useful pytest flags are:

* `-s`: disable output capture so, for example, `print` messages are shown.
* `-v`: verbose mode, displays names of test functions run.
* `-k`: run a specific test function.

Specific subsets of tests, or individual tests, can be rerun, for example:

```console
$ pytest -s -k "test_output_bam" riboviz/test/regression/test_vignette.py --expected=<DIRECTORY> --skip-workflow

$ pytest -v -k test_output_bam\[WT3AT\] riboviz/test/regression/test_vignette.py --expected=$HOME/vignette-20190802-4f28def --skip-workflow
```

Note that in such cases, if `--skip-workflow` is not provided then the workflow will be run in its entirety.

A complete run may look like the following:

```console
$ pytest -v riboviz/test/regression/test_vignette.py --expected=$HOME/vignette-20190802-4f28def --skip-workflow
============================= test session starts ==============================
...
collected 51 items                                                             

riboviz/test/regression/test_vignette.py::test_index[YAL_CDS_w_250-1] PASSED               [  1%]
riboviz/test/regression/test_vignette.py::test_index[YAL_CDS_w_250-2] PASSED               [  3%]
riboviz/test/regression/test_vignette.py::test_index[YAL_CDS_w_250-3] PASSED               [  5%]
riboviz/test/regression/test_vignette.py::test_index[YAL_CDS_w_250-4] PASSED               [  7%]
riboviz/test/regression/test_vignette.py::test_index[YAL_CDS_w_250-5] PASSED               [  9%]
riboviz/test/regression/test_vignette.py::test_index[YAL_CDS_w_250-6] PASSED               [ 11%]
riboviz/test/regression/test_vignette.py::test_index[YAL_CDS_w_250-7] PASSED               [ 13%]
riboviz/test/regression/test_vignette.py::test_index[YAL_CDS_w_250-8] PASSED               [ 15%]
riboviz/test/regression/test_vignette.py::test_index[yeast_rRNA-1] PASSED                  [ 17%]
riboviz/test/regression/test_vignette.py::test_index[yeast_rRNA-2] PASSED                  [ 19%]
riboviz/test/regression/test_vignette.py::test_index[yeast_rRNA-3] PASSED                  [ 21%]
riboviz/test/regression/test_vignette.py::test_index[yeast_rRNA-4] PASSED                  [ 23%]
riboviz/test/regression/test_vignette.py::test_index[yeast_rRNA-5] PASSED                  [ 25%]
riboviz/test/regression/test_vignette.py::test_index[yeast_rRNA-6] PASSED                  [ 27%]
riboviz/test/regression/test_vignette.py::test_index[yeast_rRNA-7] PASSED                  [ 29%]
riboviz/test/regression/test_vignette.py::test_index[yeast_rRNA-8] PASSED                  [ 31%]
riboviz/test/regression/test_vignette.py::test_tmp_fq[WT3AT-nonrRNA] PASSED                [ 33%]
riboviz/test/regression/test_vignette.py::test_tmp_fq[WT3AT-trim] PASSED                   [ 35%]
riboviz/test/regression/test_vignette.py::test_tmp_fq[WT3AT-unaligned] PASSED              [ 37%]
riboviz/test/regression/test_vignette.py::test_tmp_fq[WTnone-nonrRNA] PASSED               [ 39%]
riboviz/test/regression/test_vignette.py::test_tmp_fq[WTnone-trim] PASSED                  [ 41%]
riboviz/test/regression/test_vignette.py::test_tmp_fq[WTnone-unaligned] PASSED             [ 43%]
riboviz/test/regression/test_vignette.py::test_tmp_sam[WT3AT-orf_map_clean] PASSED         [ 45%]
riboviz/test/regression/test_vignette.py::test_tmp_sam[WT3AT-orf_map] PASSED               [ 47%]
riboviz/test/regression/test_vignette.py::test_tmp_sam[WT3AT-rRNA_map] PASSED              [ 49%]
riboviz/test/regression/test_vignette.py::test_tmp_sam[WTnone-orf_map_clean] PASSED        [ 50%]
riboviz/test/regression/test_vignette.py::test_tmp_sam[WTnone-orf_map] PASSED              [ 52%]
riboviz/test/regression/test_vignette.py::test_tmp_sam[WTnone-rRNA_map] PASSED             [ 54%]
riboviz/test/regression/test_vignette.py::test_output_bai[WT3AT] PASSED                    [ 56%]
riboviz/test/regression/test_vignette.py::test_output_bai[WTnone] PASSED                   [ 58%]
riboviz/test/regression/test_vignette.py::test_output_bam[WT3AT] PASSED                    [ 60%]
riboviz/test/regression/test_vignette.py::test_output_bam[WTnone] PASSED                   [ 62%]
riboviz/test/regression/test_vignette.py::test_output_bedgraph[WT3AT-minus] PASSED         [ 64%]
riboviz/test/regression/test_vignette.py::test_output_bedgraph[WT3AT-plus] PASSED          [ 66%]
riboviz/test/regression/test_vignette.py::test_output_bedgraph[WTnone-minus] PASSED        [ 68%]
riboviz/test/regression/test_vignette.py::test_output_bedgraph[WTnone-plus] PASSED         [ 70%]
riboviz/test/regression/test_vignette.py::test_output_h5[WT3AT] PASSED                     [ 72%]
riboviz/test/regression/test_vignette.py::test_output_h5[WTnone] PASSED                    [ 74%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WT3AT-3nt_periodicity] PASSED    [ 76%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WT3AT-codon_ribodens] PASSED     [ 78%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WT3AT-pos_sp_nt_freq] PASSED     [ 80%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WT3AT-pos_sp_rpf_norm_reads] PASSED [ 82%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WT3AT-read_lengths] PASSED       [ 84%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WT3AT-tpms] PASSED               [ 86%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WTnone-3nt_periodicity] PASSED   [ 88%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WTnone-codon_ribodens] PASSED    [ 90%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WTnone-pos_sp_nt_freq] PASSED    [ 92%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WTnone-pos_sp_rpf_norm_reads] PASSED [ 94%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WTnone-read_lengths] PASSED      [ 96%]
riboviz/test/regression/test_vignette.py::test_output_tsv[WTnone-tpms] PASSED              [ 98%]
riboviz/test/regression/test_vignette.py::test_output_tpms_collated_tsv PASSED             [100%]

========================= 51 passed in 614.86 seconds ==========================
```


---

## Run `prep_riboviz.py` tests

Run tests of error conditions and exit codes:

```console
$ pytest -v riboviz/test/tools/test_prep_riboviz.py
============================= test session starts ==============================
...

riboviz/test/tools/test_prep_riboviz.py::test_config_error_missing_config_file PASSED [ 20%]
riboviz/test/tools/test_prep_riboviz.py::test_index_error_missing_fa PASSED    [ 40%]
riboviz/test/tools/test_prep_riboviz.py::test_no_samples_error PASSED          [ 60%]
riboviz/test/tools/test_prep_riboviz.py::test_samples_error_missing_samples PASSED [ 80%]
riboviz/test/tools/test_prep_riboviz.py::test_config_error_missing_dir_in PASSED [100%]

=========================== 5 passed in 0.71 seconds ===========================
```

---

## Logging

`prep_riboviz.py` logging is handled via Python's `logging` module. The configuration file is in `riboviz/logging.yaml`.

A custom configuration file can be provided by defining a `RIBOVIZ_LOG_CONFIG` environment variable, for example:

```console
$ RIBOVIZ_LOG_CONFIG=custom_logging.yaml
```

---

## Create simulated FASTQ files

`riboviz/tools/create_fastq_examples.py` creates simple simulated FASTQ files to test adaptor trimming, UMI extraction and deduplication.

To run:

```console
$ python -m riboviz.tools.create_fastq_examples DIRECTORY
```

Or:

```console
$ python -m riboviz/tools/create_fastq_examples.py DIRECTORY
```

where `DIRECTORY` is the directory into which the simulated files are to be written. The following files are created:                                 

* `simdata_UMI5and3_4nt_adaptor.fastq`: FASTQ file with 9 reads, each with a 4nt UMI at the 5' end, a 4nt UMI at the 3' end and a 11nt adaptor at the 3' end. Reads can be grouped by UMI into 5 groups.
* `simdata_UMI5and3_4nt.fastq`: FASTQ file identical to the above but with the adaptor trimmed.
* `simdata_extracted_UMI5and3_4nt.fastq`: FASTQ file identical to the above but with the UMIs extracted and concatenated to the header.
* `simdata_extracted_UMI3_4nt.fastq`: FASTQ file with 8 reads, each with a 4nt UMI at the 3' end. Reads can be grouped by UMI into 4 groups.
* `simdata_UMI3_4nt.fastq`: FASTQ file identical to the above but with the UMI extracted and concatenated to the header.

The files with these names in `data/` were created using this script.

---

## Coding style

For Python code:

Regularly run `pylint`, `pycodestyle`, and `2to3` and update the code to adopt the recommendations as far as possible. For example, to run these on `validation.py`:

```console
$ pylint riboviz/validation.py
$ pycodestyle riboviz/validation.py
$ 2to3 riboviz/validation.py
```

For R code:

Follow [Google fork](https://google.github.io/styleguide/Rguide.html) of [Tidyverse R Style guide](https://style.tidyverse.org/) where possible - the Google fork is largely the same but differentiates more between functions (BigCamelCase) and variable names (snake_case) and does not assign to the right, for example. The Tidyverse R Style guide is implemented in `lintr` and `styleR` packages, so please bear in mind the small number of changes needed to follow the Google style.  

`Lintr` package can be used to produce programmatic output of style issues but does not edit the code, whilst the `styleR` package makes automatic adjustments to selections or files by default.

Lintr can be used within an IDE such as RStudio via an add-in once installed if preferred and run on the current file. It can also be run within the R terminal (for example on `generate_stats_figs.R`) with the command: `lint("$HOME/RiboViz/rscripts/generate_stats_figs.R")`, but if there is considerable output or you wish to work through the output bit by bit, it's possible to send it to an output file using `sink()` as below:

```R
$ R
> sink('lintR-output.txt')
> lint("$HOME/RiboViz/rscripts/generate_stats_figs.R")
> sink()
```

`StyleR` also has an add-in for RStudio IDE, which allows selected code, current file or current package to be styled automatically according to the Tidyverse style guide. It can also be run from command line (see package information for more details). There are a considerable number of options for setting 'strictness' of adherence to the style guide, which may be useful.

---

## Handling missing configuration values

### YAML `NULL` and Python `None`

If a parameter in a YAML file has value `null`, `NULL` or no value at all then, after reading the file into Python (using the `yaml` library), it will have value `None`. For example, given a YAML file with:

```yaml
a: null
b: NULL
c:
```

Loading this with

```python
with open("config.yml") as f:
    config = yaml.load(f, Loader=yaml.SafeLoader)
```

would result in the following all having value `None`:

```python
config['a']
config['b']
config['c']
```

### Check for missing or `None` values in Python dictionaries

A common pattern to check for missing or `None` values is:

```python
if 'some_key' not in config or config['some_key'] is None:
    # Some action if configuration is missing
    # e.g. set a default value, for optional configuration
    # e.g. raise an error, for mandatory configuration
```

### `NULL` vs. `NA` in R command-line parsing

When defining command-line parameters for R scripts using `make_option`, one can specify a default value of `NULL` or `NA` for optional parameters. Each of these behaves in a different way.

Consider a script which defines the following command-line options, all of which declare `default=NULL`:

```R
make_option("--t_rna", type="character", default=NULL),
make_option("--codon_pos", type="character", default=NULL),
make_option("--orf_gff_file", type="character", default=NULL),
make_option("--features_file", type="character", default=NULL),
```

and which then reads in the options into an `opt` variable and prints this variable:

```R
opt <- parse_args(OptionParser(option_list=option_list))
attach(opt)
print("Running with parameters:")
opt
```

If we run this script with one of these arguments, we see that only that argument is present in `opt`:

```console
$ Rscript script.R --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3

Running with parameters:

$orf_gff_file
[1] "vignette/input/yeast_YAL_CDS_w_250utrs.gff3"
```

Suppose, however, we declare `default=NA`:

```R
  make_option("--t_rna", type="character", default=NA),
  make_option("--codon_pos", type="character", default=NA),
  make_option("--orf_gff_file", type="character", default=NA),
  make_option("--features_file", type="character", default=NA),
```

If we now run the script, we see that all the arguments are present in `opt`:

```console
$ Rscript script.R --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3

Running with parameters:

$t_rna
[1] NA

$codon_pos
[1] NA

$orf_gff_file
[1] "vignette/input/yeast_YAL_CDS_w_250utrs.gff3"

$features_file
[1] NA
```

i.e. the options are present but with value `NA`. When defining optional command-line parameters in R, we recommend the use of `default=NA`, and not `default=NULL`, so the option is in the variable holding the parsed options.

If refactoring existing code then calls to `is.null` can be replaced by calls to `is.na`.

---

## Repository structure

```
data/            # Data files used by scripts and vignette
docs/            # Documentation
install/         # Bash scripts to install dependencies
riboviz/         # Python package code
  tools/         # End-user scripts, including prep_riboviz.py
  test/          # pytest-compliant tests
    regression/  # prep_riboviz.py and vignette regression test
rmarkdown/       # Rmarkdown scripts for data preprocessing
rscripts/        # R scripts invoked by vignette
vignette/        # Vignette configuration and input data
website/         # RiboViz Shiny server code and data
```

---

## Debugging R scripts with appropriate command-line arguments

To debug R scripts such as `generate_stats_figs.R` and `bam_to_h5.R`, they need to be run with the correct command-line arguments to discover the bug. R has good tools for interactive debugging, [explained in Hadley Wickham's chapter on debugging in R](https://adv-r.hadley.nz/debugging.html). However, interactive debugging tools such as `browser()` don't interrupt a call to `Rscript`. Instead you need to modify the call from
```
Rscript code_to_debug.R --myarg1 value1
```
to
```
R --args --myarg1 value1
```
then, from the R prompt run
```
> source('code_to_debug.R')
```
this will accept `debug()` and `browser()` statements run from the interactive R prompt.

For example, in the vignette we call:
```
Rscript --vanilla rscripts/generate_stats_figs.R --Ncores=1 --MinReadLen=10 --MaxReadLen=50 --Buffer=250 --PrimaryID=Name --dataset=vignette --hdFile=vignette/output/WTnone.h5 --out_prefix=vignette/output/WTnone --orf_fasta=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True --dir_out=vignette/output --do_pos_sp_nt_freq=True --t_rna=data/yeast_tRNAs.tsv --codon_pos=data/yeast_codon_pos_i200.RData --features_file=data/yeast_features.tsv --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 --count_threshold=64 --asite_disp_length_file=vignette/input/asite_disp_length_yeast_standard.txt
```
But to debug a new feature, instead run:
```
R --vanilla --args --Ncores=1 --MinReadLen=10 --MaxReadLen=50 --Buffer=250 --PrimaryID=Name --dataset=vignette --hdFile=vignette/output/WTnone.h5 --out_prefix=vignette/output/WTnone --orf_fasta=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True --dir_out=vignette/output --do_pos_sp_nt_freq=True --t_rna=data/yeast_tRNAs.tsv --codon_pos=data/yeast_codon_pos_i200.RData --features_file=data/yeast_features.tsv --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 --count_threshold=64 --asite_disp_length_file=vignette/input/asite_disp_length_yeast_standard.txt
```
then 
```
> source('rscripts/generate_stats_figs.R')
```
For example, to debug a specific line of code, you could add a `browser()` statement in the source first. Alternatively, you could copy and paste the parts of the code you wanted to run, as long as earlier dependencies are run first (packages, importing command args, function definitions).
