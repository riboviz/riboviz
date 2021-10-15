# Style guide

This page summarises style guidelines for the riboviz source code, documentation, parameters and files.

* ['riboviz'](#riboviz)
* [Python source code](#python-source-code)
  - [Python file names](#python-file-names)
  - [Python comments](#python-comments)
  - [Style checking Python code](#style-checking-python-code)
* [R source code](#r-source-code)
  - [R file names](#r-file-names)
  - [R comments](#r-comments)
  - [Style checking R code](#style-checking-r-code)
* [Configuration parameters](#configuration-parameters)
* [Command-line tools](#command-line-tools)
  - [Command-line parameters](#command-line-parameters)
  - [Python command-line tools](#python-command-line-tools)
  - [R command-line tools](#r-command-line-tools)
* [Input and output file names](#input-and-output-file-names)
* [Nextflow](#nextflow)
  - [Process names](#process-names)
  - [Channel names](#channel-names)
  - [Variable names](#variable-names)

---

## 'riboviz'

The sofware is called 'riboviz', one word, all lower-case, no hyphens.

Please do **not** refer to the software as RiboViz, Riboviz, Ribo-Viz, Ribo-viz or ribo-viz.

---

## Python source code

Python code should be formatted to conform as far as possible to the [PEP 8 -- Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/).

The tools [pylint](https://pylint.org/) and [pycodestyle](https://pycodestyle.pycqa.org/en/latest/) can help by highlighting deviations from the coding style. See [Style checking Python code](#style-checking-python-code) below.

### Python file names

Python file names must be in snake-case i.e., all lower case and delimited by underscores, not hyphens. For example, `sample_sheets.py`, `get_cds_codons.py`.

### Python comments

Python code should be commented using [doc-strings](https://www.python.org/dev/peps/pep-0257/). The doc-strings should include [Sphinx](https://www.sphinx-doc.org/en/master/)-compliant [reStructuredText](https://docutils.sourceforge.io/rst.html) markup. Examples are given below.

Modules should have a doc-string descring the module. For example:

```python
"""
Workflow configuration parameter names.
"""
```

Constants and module-level variable should have doc-strings briefly describing what they are. For example:

```python
INPUT_DIR = "dir_in"
""" Input directory. """
```

Functions should have doc-strings describing what they do and their parameters, return values (if applicable) and any exceptions (if applicable). Types of parameters and return values should also be documented. For example:

```python
def load_sample_sheet(file_name, delimiter="\t", comment="#"):
    """
    Load a sample sheet from a file. The sample sheet is assumed to
    have a header with column names ``SampleID`` and ``TagRead``.

    :param file_name: File name
    :type file_name: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    :param comment: Comment prefix
    :type comment: str or unicode
    :return: Sample sheet
    :rtype: pandas.core.frame.DataFrame
    :raise FileNotFoundError: If the file cannot be found or is \
    not a file
    :raise AssertionError: If there is no header with ``SampleID`` \
    and ``TagRead`` columns.
    """
```

To use fixed-width (teletype-style) text use ` `` ` (see the markup on `SampleID` and `TagRead` above).

For references to other Python modules, functions or constants use Sphinx's `:py:` roles. For example:

* Cross-reference to a module (`:py:mod:`):

```python
Count reads using :py:mod:`riboviz.tools.count_reads`.
```

* Cross-reference to a constant (`:py:constant:`):

```python
... also matching :py:const:`riboviz.workflow_files.UMI_EXTRACT_FQ` ...
```

* Cross-reference to a function (`:py:func:`):

```python
See :py:func:`riboviz.count_reads.count_reads` for information ...
```

Note that if the cross-reference is to an entity in the same module then only the local name of the entity needs to be specified. For example:

```python
See :py:func:`count_reads` for information ...
```

For examples and further information, see:

* Sphinx's Python [Info field lists](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#info-field-lists) on the Python fields that Sphinx can recognise.
* Sphinx's [Cross-referencing Python objects](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#cross-referencing-python-objects) on the available cross-references.
* Sphinx's [The Python Domain](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#the-python-domain).

### Style checking Python code

[pylint](https://pylint.org/) and [pycodestyle](https://pycodestyle.pycqa.org/en/latest/) can be run on individual files, groups of files or every file in a directory and its subdirectories. For example:

```console
$ pylint riboviz/check_fasta_gff.py
$ pycodestyle riboviz/check_fasta_gff.py
```
```console
$ pylint riboviz/
$ pycodestyle riboviz/
```

---

## R source code

R code should be commented using [Google’s R Style Guide](https://google.github.io/styleguide/Rguide.html), a fork of [Tidyverse R Style guide](https://style.tidyverse.org/), where possible. The Google fork is largely the same as Tidyverse but differentiates more between functions (`BigCamelCase`) and variable names (`snake_case`) and does not assign to the right, for example.

The packages [lintr](https://cran.r-project.org/package=lintr) and [styler](https://styler.r-lib.org/) can help by highlighting deviations from the coding style. See [Style checking R code](#style-checking-r-code) below.

### R file names

R file names must be in snake-case i.e., all lower case and delimited by underscores, not hyphens. For example, `generate_stats_figs.R`.

### R Comments

R code should be commented using [roxygen2](https://roxygen2.r-lib.org/)-compliant comments.

See roxygen2's [Rd formatting](https://roxygen2.r-lib.org/articles/rd-formatting.html) for a comprehensive guide.

### Style checking R code

[lintr](https://cran.r-project.org/package=lintr) can be used to produce a programmatic output of style issues.

lintr can be used within an IDE such as RStudio via an add-in, once installed, if preferred, and run on the current file. lintr can also be run within the R terminal. For example, it can be run on `generate_stats_figs.R` with the command:
```console
$ R
```
```R
> library(lintr)
> lint(rscripts/generate_stats_figs.R")
```

If there is considerable output or you wish to work through the output bit by bit, it is possible to save it to an output file, using `sink()`, as follows:

```console
$ R
```
```R
> sink('lintR-output.txt')
> lint("$HOME/riboviz/rscripts/generate_stats_figs.R")
> sink()
```

[styler](https://styler.r-lib.org/) also has an add-in for RStudio, which allows the selected code, current file or current package to be styled automatically according to the Tidyverse style guide. It can also be run from command line (see package information for more details). There are a considerable number of options for setting 'strictness' of adherence to the style guide, which may be useful.

**Note:**

As commented above, `lint` checks code for conformance to the Tidyverse R style guide which recommends `snake_case` function names. However, we prefer Google-style `CamelCase` function names. This means you can expect to see the following types of warnings when running `lint` on riboviz R code:

```
rscripts/generate_stats_figs.R:179:1: style: Variable and function name style should be snake-case.
ThreeNucleotidePeriodicity <- function(gene_names, dataset, hd_file, gff_df) {
^~~~~~~~~~~~~~~~~~~~~~~~~~
```

---

## Configuration parameters

Configuration parameters must be in snake-case i.e., all lower case and delimited by underscores, not hyphens. For example, `adapters` or `asite_disp_length_file` or `orf_gff_file`.

---

## Command-line tools

### Command-line parameters

For consistency with bash and other command-line tools, command-line parameters should be implemented as one, or both, of:

* A single alphanumeric character prefixed by a hyphen. For example, `-v`, `-c 123`, `-s GATC`.
* A kebab-case token i.e., all lower case and delimited by hyphens, not underscores, and prefixed by two hyphens. For example, `--verbose`, `--control=123`, `--match-sequence=GATC`.

### Python command-line tools

Python command-line tools must be placed within `riboviz/tools/` i.e. within the `riboviz.tools` module.

Command-line parsing must be implemented using [argparse](https://docs.python.org/3/library/argparse.html).

**Defining command-line parameters**

The `dest` parameter of `ArgumentParser.add_argument` can be use to explicitly define the Python variables into which a command-line parameter is to be placed. For example:

```python
parser = argparse.ArgumentParser(description="Program")
parser.add_argument("-o", "--output-dir", dest="output_dir", nargs='?',
                    help="Output directory")
options = parser.parse_args()
output_dir = options.output_dir
```

### R command-line tools

R command-line tools must be placed within `rscripts/`.

Command-line parsing must be implemented using [optparse](https://cran.r-project.org/web/packages/optparse/index.html).

**Defining command-line parameters**

You can call `parse_args` with `convert_hyphens_to_underscores=TRUE` to automatically convert hyphens in command-line parameters to underscores. For example, the following shows how an option with a hyphen is accessed within R as a variable with an underscore:

```R
option_list <- list( 
  make_option("--output-dir", type="character", default="./",
              help="Output directory"))
parser <- OptionParser(option_list=option_list)
opts <- parse_args(parser,
                   convert_hyphens_to_underscores=TRUE)
output_dir <- opts$options$output_dir
```

**Use `default=NA` not `default=NULL`**

If defining optional parameters that take no value by default then use `default=NA`, and **not** `default=NULL`, so that the option is in the variable holding the parsed options.

If this is not done then the option will not be present in the variable holding the parsed options. For example, consider a script which defines the following command-line options, all of which declare `default=NULL`:

```R
make_option("--t-rna-file", type="character", default=NULL),
make_option("--codon-positions-file", type="character", default=NULL),
make_option("--orf-gff-file", type="character", default=NULL),
make_option("--features-file", type="character", default=NULL),
```

If this script then reads in the options into an `opt` variable and prints this variable:

```R
opt <- parse_args(OptionParser(option_list=option_list),
                  convert_hyphens_to_underscores=TRUE)
attach(opt)
print("Running with parameters:")
opt
```

then, if we run this script with only one of these arguments, we see that only that argument is present in `opt`:

```console
$ Rscript script.R --orf-gff-file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3

Running with parameters:

$orf_gff_file
[1] "vignette/input/yeast_YAL_CDS_w_250utrs.gff3"
```

If, however, we declare `default=NA`, for each parameter:

```R
make_option("--t-rna-file", type="character", default=NA),
make_option("--codon-positions-file", type="character", default=NA),
make_option("--orf-gff-file", type="character", default=NA),
make_option("--features-file", type="character", default=NA),
```

then, if we now run the script, we see that all the arguments are present in `opt`:

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

i.e. the options are present but with value `NA`.

The function `is.na` (as opposed to `is.null`) can be used to check if a variable has value `NA`.

---

## Input and output file names

riboviz-specific input and output file names must be in snake-case i.e., lower case and delimited by underscores, not hyphens. Upper-case is permitted for acronyms e.g., `ORF`, `CDS`, `APE`, `TPMs`, `RNA`.

For example, `vignette_config.yaml`, `read_counts_per_file.tsv`, `TPMs_all_CDS_all_samples.tsv`.

---

## Nextflow

### Process names

Nextflow process names must be in CamelCase with the first letter being lower-case. Upper-case is peromitted for acronyms e.g., `ORF`, `CDS`, `APE`, `TPMs`, `RNA`.

For example, `buildIndicesrRNA`, `demultiplex`, `staticHTML`.

### Channel names

Nextflow channel names must be in snake-case i.e., lowe case and delimited by underscores, not hyphens.

For channels which are file names ensure the channel name includes the file type as its last component, delimited by an underscore.

For example, `multiplex_sample_sheet_tsv`, `sample_fq`

### Variable names

Nextflow variable names must be in snake-case i.e., lowe case and delimited by underscores, not hyphens.
