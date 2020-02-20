# Coding style

## Python

Regularly run `pylint`, `pycodestyle`, and `2to3` and update the code to adopt the recommendations as far as possible. For example, to run these on `riboviz/compare_files.py`:

```console
$ pylint riboviz/compare_files.py
$ pycodestyle riboviz/compare_files.py
$ 2to3 riboviz/compare_files.py
```

---

## R

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
