# Developing R components

* [R style](#r-style)
  - [File names](#file-names)
  - [Comments](#comments)
  - [Style checking R code](#style-checking-r-code)
* [R command-line tools](#r-command-line-tools)
  - [Defining command-line parameters](#defining-command-line-parameters)
  - [Use `default=NA` not `default=NULL`](#use-defaultna-not-defaultnull)
* [Debugging R scripts with appropriate command-line arguments](#debugging-r-scripts-with-appropriate-command-line-arguments)
* [Run R tests](#run-r-tests)

---

## R style

R code should be commented using [Google’s R Style Guide](https://google.github.io/styleguide/Rguide.html), a fork of [Tidyverse R Style guide](https://style.tidyverse.org/), where possible. The Google fork is largely the same as Tidyverse but differentiates more between functions (`BigCamelCase`) and variable names (`snake_case`) and does not assign to the right, for example.

The packages [lintr](https://cran.r-project.org/package=lintr) and [styler](https://styler.r-lib.org/) can help by highlighting deviations from the coding style. See [Style checking R code](#style-checking-r-code) below.

See also the riboviz [Style guide](./style-guide.md).

### File names

R file names must be in snake-case i.e., lower-case and delimited by underscores, not hyphens. For example, `generate_stats_figs.R`.

### Comments

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

## R command-line tools

R command-line tools must be placed within `rscripts/`.

Command-line parsing must be implemented using [optparse](https://cran.r-project.org/web/packages/optparse/index.html).

### Defining command-line parameters

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

### Use `default=NA` not `default=NULL`

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

## Debugging R scripts with appropriate command-line arguments

To debug R scripts such as `generate_stats_figs.R` and `bam_to_h5.R`, they need to be run with the correct command-line arguments to discover the bug. R has good tools for interactive debugging, [explained in Hadley Wickham's chapter on debugging in R](https://adv-r.hadley.nz/debugging.html). However, interactive debugging tools such as `browser()` don't interrupt a call to `Rscript`. Instead you need to modify the call from:

```console
$ Rscript code_to_debug.R --myarg1 value1
```

to:

```console
$ R --args --myarg1 value1
```

then, from the R prompt run:

```R
> source('code_to_debug.R')
```

this will accept `debug()` and `browser()` statements run from the interactive R prompt.

For example, in the vignette we call:

```console
$ Rscript --vanilla /home/ubuntu/riboviz/rscripts/generate_stats_figs.R \
  --num-processes=1            --min-read-length=10 \
  --max-read-length=50            --buffer=250 \
  --primary-id=Name            --dataset=vignette \
  --hd-file=WTnone.h5 \
  --orf-fasta-file=yeast_YAL_CDS_w_250utrs.fa            --rpf=true \
  --output-dir=.            --do-pos-sp-nt-freq=true \
  --t-rna-file=yeast_tRNAs.tsv \
  --codon-positions-file=yeast_codon_pos_i200.RData \
  --features-file=yeast_features.tsv \
  --orf-gff-file=yeast_YAL_CDS_w_250utrs.gff3 \
  --asite-disp-length-file=yeast_standard_asite_disp_length.txt \
  --count-threshold=64
```

But to interactively debug a new feature, we'd run:

```console
$ R --vanilla --args \
  --num-processes=1            --min-read-length=10 \
  --max-read-length=50            --buffer=250 \
  --primary-id=Name            --dataset=vignette \
  --hd-file=WTnone.h5 \
  --orf-fasta-file=yeast_YAL_CDS_w_250utrs.fa            --rpf=true \
  --output-dir=.            --do-pos-sp-nt-freq=true \
  --t-rna-file=yeast_tRNAs.tsv \
  --codon-positions-file=yeast_codon_pos_i200.RData \
  --features-file=yeast_features.tsv \
  --orf-gff-file=yeast_YAL_CDS_w_250utrs.gff3 \
  --asite-disp-length-file=yeast_standard_asite_disp_length.txt \
  --count-threshold=64
```

then:

```R
> source('<PATH_TO_RIBOVIZ_DIRECTORY>/rscripts/generate_stats_figs.R')
```

To debug a specific line of code, you could add a `browser()` statement in the source first. Alternatively, you could copy and paste the parts of the code you wanted to run, as long as earlier dependencies are run first (packages, importing command arguments, function definitions).

**Note:** at present, the riboviz R scripts `bam_to_h5.R`, `generate_stats_figs.R` and `collate_tpms.R` import other riboviz R scripts. If running these riboviz R scripts interactively, via `R` and `source`, then the directory in which they are run must be such that `rscripts` is a sibling of an ancestor of the directory in which the script is run interactively. For example, running a script interactively within a sub-sub-directory of Nextflow's `work/` directory or a `debug_gen_stats_figs` directory (as described in the next section) where either of these directories are in the same directory as `rscripts`.

---

## Run R tests

testthat-compliant tests for R scripts can be run as follows:

```console
$ Rscript rscripts/tests/testthat.R
```

Individual test scripts can be run as follows (for example):

```console
$ Rscript rscripts/tests/testthat/test_bam_to_h5.R 
```

Test scripts can also be run from within R as follows (for example):

```R
> library(testthat)
> test_file("rscripts/tests/testthat/test_bam_to_h5.R")
```
