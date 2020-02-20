# Coding style

## Documenting code

### Python

Python code should be documented using [doc-strings](https://www.python.org/dev/peps/pep-0257/). The doc-strings should include [Sphinx](https://www.sphinx-doc.org/en/master/)-compliant [reStructuredText](https://docutils.sourceforge.io/rst.html) markup. Examples are given below.

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

Functions should have doc-strings describing what they do and their parameters, return values (if applicable) and any exceptions (if applicable). Types of parameters and return values should also be documented. These should be documented using specific fields. For example:

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

Note that if the cross-reference is to an entity in the same module only the local name of the entity needs to be specified. For example:

```python
See :py:func:`count_reads` for information ...
```

For examples and further information, see:

* Our source code.
* Sphinx's Python [Info field lists](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#info-field-lists) on the Python fields that Sphinx can recognise.
* Sphinx's [Cross-referencing Python objects](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#cross-referencing-python-objects) on the available cross-references.
* Sphinx's [The Python Domain](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#the-python-domain)

---

## Style checking

### Python

Regularly run `pylint`, `pycodestyle`, and `2to3` and update the code to adopt the recommendations as far as possible. For example, to run these on `riboviz/compare_files.py`:

```console
$ pylint riboviz/compare_files.py
$ pycodestyle riboviz/compare_files.py
$ 2to3 riboviz/compare_files.py
```

### R

Follow [Google’s R Style Guide](https://google.github.io/styleguide/Rguide.html), a fork of [Tidyverse R Style guide](https://style.tidyverse.org/), where possible - the Google fork is largely the same but differentiates more between functions (BigCamelCase) and variable names (snake_case) and does not assign to the right, for example. The Tidyverse R Style guide is implemented in `lintr` and `styleR` packages, so please bear in mind the small number of changes needed to follow the Google style.  

`Lintr` package can be used to produce programmatic output of style issues but does not edit the code, whilst the `styleR` package makes automatic adjustments to selections or files by default.

Lintr can be used within an IDE such as RStudio via an add-in once installed if preferred and run on the current file. It can also be run within the R terminal. For example, it can be run on `generate_stats_figs.R` with the command:

```console
$ R
```
```R
> lint("$HOME/RiboViz/rscripts/generate_stats_figs.R")
```

If there is considerable output or you wish to work through the output bit by bit, it's possible to send it to an output file using `sink()` as below:

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

### Python

If using Python's [argparse](https://docs.python.org/3/library/argparse.html) package, the `dest` parameter of `ArgumentParser.add_argument` can be use to explicitly define the Python variables into which a command-line is to be placed. For example:

```python
parser = argparse.ArgumentParser(description="Program")
parser.add_argument("-o", "--output-dir", dest="output_dir", nargs='?',
                    help="Output directory")
options = parser.parse_args()
output_dir = options.output_dir
```

### R

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
