# Developer guide

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

`riboviz.tools.compare_files` is a script that can compare the files output by any stage of the workflow, `riboviz.tools.prep_riboviz`.

It can be run as follows:

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

---

## Run vignette regression test suite

The vignette regression test suite optionally runs `prep_riboviz` on the vignette data, in `vignette/`, then compares the results, in `vignette/`, to a directory of pre-calculated results.

The tests can be run using pytest:

```console
$ pytest riboviz/test/regression/test_vignette.py \
        --expected=<DIRECTORY> \
        [--skip-workflow] \
	[--check-index-tmp]
```

where:

* `--expected=<DIRECTORY>`: directory with expected vignette files. This is assumed to have `index/` `tmp/` and `output/` directories.
* `--skip-workflow`: request that the `prep_riboviz` workflow **not** be run, instead use existing data files in `vignette/` for testing.
* `--check-index-tmp`: request that the index and temporary files be checked also (the default is that only the output files are checked)

For each file output by the vignette, the files are compared using the same comparisons as for `compare_files` above.

If `--check-index-tmp` is not provided then tests for index and temporary files will be skipped. This will appear as follows:

```console
riboviz/test/regression/test_vignette.py ssssssssssssssssssssssssssssss.....
```

or, if using pytest's `-v`, verbose, mode:

```console
...
riboviz/test/regression/test_vignette.py::test_index[1-YAL_CDS_w_250] SKIPPED [  1%]
 riboviz/test/regression/test_vignette.py::test_index[1-yeast_rRNA] SKIPPED [  2%]
riboviz/test/regression/test_vignette.py::test_index[2-YAL_CDS_w_250] SKIPPED [  4%]
riboviz/test/regression/test_vignette.py::test_index[2-yeast_rRNA] SKIPPED [  5%]
...
```

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

---

## Run tests

Run all tests (excluding regression tests):

```console
$ pytest -v --ignore-glob="*regression*"
```

Run all tests and generate test coverage report (including lines that were not invoked):

```console
$ pytest --cov-config=.coveragerc --cov-report term-missing \
    --cov=riboviz --ignore-glob="*regression*"
```

---

## Logging

`prep_riboviz` logging is handled via Python's `logging` module. The configuration file is in `riboviz/logging.yaml`.

A custom configuration file can be provided by defining a `RIBOVIZ_LOG_CONFIG` environment variable, for example:

```console
$ RIBOVIZ_LOG_CONFIG=custom_logging.yaml
```

---

## Coding style

For Python code:

Regularly run `pylint`, `pycodestyle`, and `2to3` and update the code to adopt the recommendations as far as possible. For example, to run these on `riboviz/compare_files.py`:

```console
$ pylint riboviz/compare_files.py
$ pycodestyle riboviz/compare_files.py
$ 2to3 riboviz/compare_files.py
```

For R code:

Follow [Google fork](https://google.github.io/styleguide/Rguide.html) of [Tidyverse R Style guide](https://style.tidyverse.org/) where possible - the Google fork is largely the same but differentiates more between functions (BigCamelCase) and variable names (snake_case) and does not assign to the right, for example. The Tidyverse R Style guide is implemented in `lintr` and `styleR` packages, so please bear in mind the small number of changes needed to follow the Google style.  

`Lintr` package can be used to produce programmatic output of style issues but does not edit the code, whilst the `styleR` package makes automatic adjustments to selections or files by default.

Lintr can be used within an IDE such as RStudio via an add-in once installed if preferred and run on the current file. It can also be run within the R terminal (for example on `generate_stats_figs.R`) with the command: `lint("$HOME/RiboViz/rscripts/generate_stats_figs.R")`, but if there is considerable output or you wish to work through the output bit by bit, it's possible to send it to an output file using `sink()` as below:

```console
$ R
```
```R
> sink('lintR-output.txt')
> lint("$HOME/RiboViz/rscripts/generate_stats_figs.R")
> sink()
```

`StyleR` also has an add-in for RStudio IDE, which allows selected code, current file or current package to be styled automatically according to the Tidyverse style guide. It can also be run from command line (see package information for more details). There are a considerable number of options for setting 'strictness' of adherence to the style guide, which may be useful.

---

## Command-line parameters

For consistency with bash and other command-line tools, command-line parameters should be implemented one, or both, of:

* A single alphanumeric character prefixed by a hyphen e.g. `-v`, `-c 123`, `-s GATC`.
* A sequence of lower-case alphanumeric characters, prefixed by two hyphens and delimited by hyphens e.g. `--verbose`, `--control=123`, `--match-sequence=GATC`.

If using Python's [argparse](https://docs.python.org/3/library/argparse.html) package, the `dest` parameter of `ArgumentParser.add_argument` can be use to explicitly define the Python variables into which a command-line is to be placed. For example:

```python
parser = argparse.ArgumentParser(description="Program")
parser.add_argument("-o", "--output-dir", dest="output_dir", nargs='?',
                    help="Output directory")
options = parser.parse_args()
output_dir = options.output_dir
```    

If using R's [optparse](https://cran.r-project.org/web/packages/optparse/index.html) package, call `parse_args` with `convert_hyphens_to_underscores=TRUE` to automatically convert hyphens in command-line parameters to underscores. For example, the following shows how an option with a hyphen is accessed within R as a variable with an underscore:

```R
option_list <- list( 
  make_option("--output-dir", type="character", default="./",
              help="Output directory"))
parser <- OptionParser(option_list=option_list)
opts <- parse_args(parser,
                   convert_hyphens_to_underscores=TRUE)

output_dir <- opts$options$output_dir
```

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
make_option("--t-rna-file", type="character", default=NULL),
make_option("--codon-positions-file", type="character", default=NULL),
make_option("--orf-gff-file", type="character", default=NULL),
make_option("--features-file", type="character", default=NULL),
```

and which then reads in the options into an `opt` variable and prints this variable:

```R
opt <- parse_args(OptionParser(option_list=option_list),
                  convert_hyphens_to_underscores=TRUE)
attach(opt)
print("Running with parameters:")
opt
```

If we run this script with one of these arguments, we see that only that argument is present in `opt`:

```console
$ Rscript script.R --orf-gff-file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3

Running with parameters:

$orf_gff_file
[1] "vignette/input/yeast_YAL_CDS_w_250utrs.gff3"
```

Suppose, however, we declare `default=NA`:

```R
make_option("--t-rna-file", type="character", default=NA),
make_option("--codon-positions-file", type="character", default=NA),
make_option("--orf-gff-file", type="character", default=NA),
make_option("--features-file", type="character", default=NA),
```

If we now run the script, we see that all the arguments are present in `opt`:

```console
$ Rscript script.R --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3

Running with parameters:

$t_rna_file
[1] NA

$codon_positions_file
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

```console
Rscript code_to_debug.R --myarg1 value1
```

to

```console
R --args --myarg1 value1
```

then, from the R prompt run

```R
> source('code_to_debug.R')
```

this will accept `debug()` and `browser()` statements run from the interactive R prompt.

For example, in the vignette we call:

```console
Rscript --vanilla rscripts/generate_stats_figs.R --Ncores=1 --MinReadLen=10 --MaxReadLen=50 --Buffer=250 --PrimaryID=Name --dataset=vignette --hdFile=vignette/output/WTnone.h5 --out_prefix=vignette/output/WTnone --orf_fasta=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True --dir_out=vignette/output --do_pos_sp_nt_freq=True --t_rna=data/yeast_tRNAs.tsv --codon_pos=data/yeast_codon_pos_i200.RData --features_file=data/yeast_features.tsv --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 --count_threshold=64 --asite_disp_length_file=vignette/input/asite_disp_length_yeast_standard.txt
```

But to debug a new feature, instead run:

```console
R --vanilla --args --Ncores=1 --MinReadLen=10 --MaxReadLen=50 --Buffer=250 --PrimaryID=Name --dataset=vignette --hdFile=vignette/output/WTnone.h5 --out_prefix=vignette/output/WTnone --orf_fasta=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True --dir_out=vignette/output --do_pos_sp_nt_freq=True --t_rna=data/yeast_tRNAs.tsv --codon_pos=data/yeast_codon_pos_i200.RData --features_file=data/yeast_features.tsv --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 --count_threshold=64 --asite_disp_length_file=vignette/input/asite_disp_length_yeast_standard.txt
```

then 

```R
> source('rscripts/generate_stats_figs.R')
```

For example, to debug a specific line of code, you could add a `browser()` statement in the source first. Alternatively, you could copy and paste the parts of the code you wanted to run, as long as earlier dependencies are run first (packages, importing command args, function definitions).

---

## Editing workflow images

Workflow images in `images/` are written in the open source [dot](https://graphviz.org/doc/info/lang.html) language from [GraphViz](https://www.graphviz.org/). For an overview, see [Drawing graphs with dot](https://www.graphviz.org/pdf/dotguide.pdf). GraphViz includes a command-line tool, `dot`, for converting dot files into image in various formats.

Convert `.dot` file to `.svg` file:

```console
$ dot -Tsvg workflow.dot > workflow.svg
```

When `.dot` files are updated, the corresponding `.svg` images should be updated also.

---

## Sphinx code documentation

### Updating HTML documentation

Create Sphinx pages to reference source code:

```console
$ sphinx-apidoc -f -o code-docs/ riboviz
```

Create HTML documentation

```console
$ cd code-docs
$ make html
```

Open `code-docs/_build/html/index.html` in a browser.

## How the Sphinx documentation files were originally created

```console
$ sphinx-quickstart code-docs
> Separate source and build directories (y/n) [n]: y
> Project name: RiboViz
> Author name(s): The University of Edinburgh; Rutgers University; University of California, Berkeley
> Project release []: 
> Project language [en]: 
```

Edit `code-docs/conf.py`:

* Uncomment:

```python
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
```

* Update above to:

```python
sys.path.insert(0, os.path.abspath('../../riboviz'))
```

* Update:

```python
extensions = [
]
```

* to:

```python
extensions = [
    'sphinx.ext.autodoc'
]
```

Edit `code-docs/index.rst` and replace content with:

```
RiboViz code documentation
==========================

.. toctree::
   :maxdepth: 1
   :caption: Code documentation

   modules

Indices and tables:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```
