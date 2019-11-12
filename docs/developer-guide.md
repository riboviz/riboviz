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

* `example_umi5_umi3_umi_adaptor.fastq`: FASTQ file with 9 reads, each with a 4nt UMI at the 5' end, a 4nt UMI at the 3' end and a 11nt adaptor at the 3' end. Reads can be grouped by UMI into 5 groups.
* `example_umi5_umi3_umi.fastq`: FASTQ file identical to the above but with the adaptor trimmed.                                                     
* `example_umi5_umi3.fastq`: FASTQ file identical to the above but with the UMIs extracted and concatenated to the header, with a "_" delimiter.
* `example_umi3_umi_adaptor.fastq`: FASTQ file with 8 reads, each with a 4nt UMI at the 3' end and a 11nt adaptor at the 3' end. Reads can be grouped by UMI into 4 groups.
* `example_umi3_umi.fastq`: FASTQ file identical to the above but with the adaptor trimmed.
* `example_umi3.fastq`: FASTQ file identical to the above but with the UMI extracted and concatenated to the header, with a "_" delimiter.
* `example_multiplex_umi_barcode_adaptor.fastq`: FASTQ file with 90 reads:
  - Each read has a 4nt UMI at the 5' end, a 4nt UMI at the 3' end, a 3nt barcode at the 3' end and a 11nt adaptor at the 3' end.
  - There are 9 reads for each of the following barcodes:
    - `ACG`, `ACT`, `TAG`
    - `GAC`, `GTC`, `GTA`
    - `CGA`, `TGA`, `CTT`
  - The second and third barcodes in each list have a mismatch of 1nt and 2nt respectively with the first barcode in each list.
  - When the file is demultiplexed, assuming up to 2 mismatches are allowed, then 3 sets of 27 reads will be produced, grouped by the 1st barcode in each list.
  - There are 9 reads with barcode `TTT`, which has a mismatch of 3nts to `ACG`, `GAC`, `CGA`. When the file is demultiplexed, assuming up to 2 mismatches are allowed, then these 9 reads will be unassigned.
* `example_multiplex_umi_barcode.fastq`: FASTQ file identical to the above but with the adaptor trimmed.
* `example_multiplex.fastq`: FASTQ file identical to the above but  with the barcode and UMIs extracted into the header and delimited by "_".
* `example_multiplex_tag0|1|2.fastq`: FASTQ files each with 27 reads representing the above file, demultiplexed according to the barcodes `ACG`, `GAC`, `CGA`.
* `example_multiplex_barcodes.tsv`: tab-separated values file with  `SampleID` column (with values `Tag0|1|2`) and `TagRead` column (with values `ACG`, `GAC`, `CGA`)

The files with these names in `data/example/` were created using this script.

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
#! python
"""
Creates simple simulated FASTQ files to test UMI/deduplication,
adaptor trimming, and demultiplexing.

Usage:

    python -m riboviz.tools.create_fastq_examples DIRECTORY

where `DIRECTORY` is the directory into	which the simulated files are
to be written. The following files are created:

* `example_umi5_umi3_umi_adaptor.fastq`: FASTQ file with 9 reads,
  each with a 4nt UMI at the 5' end, a 4nt UMI at the 3' end and a
  11nt adaptor at the 3' end. Reads can be grouped by UMI into 5
  groups.
* `example_umi5_umi3_umi.fastq`: FASTQ file identical to the above but
  with the adaptor trimmed.
* `example_umi5_umi3.fastq`: FASTQ file identical to the
  above but with the UMIs extracted and concatenated to the header,
  with a "_" delimiter.
* `example_umi3_umi_adaptor.fastq`: FASTQ file with 8 reads, each
  with a 4nt UMI at the 3' end and a 11nt adaptor at the 3' end. Reads
  can be grouped by UMI into 4 groups.
* `example_umi3_umi.fastq`: FASTQ file identical to the above but
  with the adaptor trimmed.
* `example_umi3.fastq`: FASTQ file identical to the above but with
   the UMI extracted and concatenated to the header, with a "_"
   delimiter.
"""

import csv
import os
import os.path
from random import choices, seed
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from riboviz.utils import BARCODE_DELIMITER
from riboviz.utils import UMI_DELIMITER


QUALITY_MEDIUM = list(range(30, 41))
""" List of medium quality scores. """
QUALITY_HIGH = list(range(39, 41))
""" List of high quality scores. """
FASTQ_FORMAT = "fastq"
""" Format string for use with Bio.SeqIO.write. """


def simulate_quality(k, qualities=QUALITY_MEDIUM, weights=None):
    """
    Simulate quality scores. This is a thin wrapper around
    random.choices whose default values represent medium Phred quality
    values.

    See https://docs.python.org/3/library/random.html#random.choices

    :param k: Number of quality scores requested
    :type k: int
    :param qualities: Available quality scores
    :type qualities: list(int)
    :param weights: Optional weightings for qualities
    :type weights: list(ints)
    :return k quality scores, selected at random from qualities
    :rtype: list(int)
    """
    return choices(qualities, k=k, weights=weights)


def make_fastq_record(name, reads, scores=None, qualities=QUALITY_MEDIUM):
    """
    Make a fastq record with sequence readall and name.

    :param name: Name
    :type name: str or unicode
    :param reads: Reads
    :type reads: str or unicode
    :param scores: Quality scores
    :type scores: list(int)
    :param qualities: Available quality scores, used if scores is
    None, in conjunction with simulate_quality to calculate quality
    scores
    :type qualities: list(int)
    :return: fastq record
    :rtype: Bio.SeqRecord.SeqRecord
    """
    if scores is None:
        scores = simulate_quality(len(reads), qualities=qualities)
    record = SeqRecord(Seq(reads, IUPAC.ambiguous_dna),
                       id=name,
                       name=name,
                       description=name)
    record.letter_annotations["phred_quality"] = scores
    return record


def trim_fastq_record_3prime(record,
                             trim,
                             add_trim=False,
                             delimiter=UMI_DELIMITER):
    """
    Copy fastq record, but trim sequence and quality scores at 3' end
    by given length.

    :param record: fastq record
    :type record: Bio.SeqRecord.SeqRecord
    :param trim: Number of nts to trim by
    :type trim: int
    :param delimiter: Delimiter to use, if add_trim is True
    :type delimiter: str or unicode
    :param add_trim: Add subsequence that was removed to record ID
    :type add_trim: bool
    :return: fastq record
    :rtype: Bio.SeqRecord.SeqRecord
    """
    quality = record.letter_annotations["phred_quality"]
    sequence = str(record.seq)
    record_extension = ""
    if add_trim:
        record_extension = delimiter + sequence[-trim:]
    return make_fastq_record(record.id + record_extension,
                             sequence[0:-trim],
                             quality[0:-trim])


def trim_fastq_record_5prime(record,
                             trim,
                             add_trim=False,
                             delimiter=UMI_DELIMITER):
    """
    Copy fastq record, but trim sequence and quality scores at 5' end
    by given length.

    :param record: fastq record
    :type record: Bio.SeqRecord.SeqRecord
    :param trim: Number of nts to trim by
    :type trim: int
    :param add_trim: Add subsequence that was removed to record ID
    :type add_trim: bool
    :param delimiter: Delimiter to use, if add_trim is True
    :type delimiter: str or unicode
    :return: fastq record
    :rtype: Bio.SeqRecord.SeqRecord
    """
    quality = record.letter_annotations["phred_quality"]
    sequence = str(record.seq)
    record_extension = ""
    if add_trim:
        record_extension = delimiter + sequence[0:trim]
    return make_fastq_record(record.id + record_extension,
                             sequence[trim:],
                             quality[trim:])


def make_fastq_records(tag,
                       read,
                       qualities,
                       umi5="",
                       umi3="",
                       barcode="",
                       adaptor="",
                       post_adaptor_nt=""):
    """
    Create a set of complementary fastq records.

    - A record for sequence: umi5 + read + umi3 + barcode + adaptor +
      post_adaptor_nt.
    - A record as above with the adaptor and post_adaptor_nt trimmed:
      umi5 + read + umi3 + barcode.
      - If adaptor and post_adaptor_nt are both "" then this is
        equivalent to the above record.
    - A record as above with the barcode and UMIs trimmed and added to
      the header with "_" delimiters: read.
      - If both the barcode and both UMIs are "" then this is
        equivalent to the above record.
      - Depending on whether barcode, umi5 and umi3 are "" the header
        will be extended with one of:
        - "<barcode>_<umi5><umi3>"
        - "<barcode>_<umi3>"
        - "_<umi5><umi3>"
        - "<umi3>"

    :param tag: Human-readable tag.
    :type tag: str or unicode
    :param read: Read
    :type read: str or unicode
    :param qualities: Available quality scores
    :type qualities: list(int)
    :param umi5: 5' end UMI
    :type umi5: str or unicode
    :param umi3: 3' end UMI
    :type umi3: str or unicode
    :param barcode: 3' end barcode
    :type barcode: str or unicode
    :param adaptor: 3' end adaptor
    :type adaptor: str or unicode
    :param post_adaptor_nt: 3' end post-adaptor nts
    :type post_adaptor_nt: str or unicode
    :returnL full record, adaptor-trimmed record, barcode- and
    UMI-extracted record
    :rtype: tuple(Bio.SeqRecord.SeqRecord, Bio.SeqRecord.SeqRecord,
    Bio.SeqRecord.SeqRecord)
    """
    sequence = umi5 + read + umi3 + barcode + adaptor + post_adaptor_nt
    record = make_fastq_record(tag, sequence, qualities=qualities)
    # Record after adaptor trimming.
    trim_record = trim_fastq_record_3prime(
        record,
        len(adaptor) + len(post_adaptor_nt))
    if barcode != "":
        # Add barcode to record ID, using "_" delimiter for consistency
        # with UMI-tools.
        barcode_ext_record = trim_fastq_record_3prime(
            trim_record,
            len(barcode),
            True,
            BARCODE_DELIMITER)
    else:
        barcode_ext_record = trim_record
    umi3_delimiter = UMI_DELIMITER
    if umi5 != "":
        # Record after 5' UMI extraction.
        # Add UMI to record ID, using "_" delimiter for consistency
        # with UMI-tools.
        umi5_ext_record = trim_fastq_record_5prime(
            barcode_ext_record,
            len(umi5),
            True,
            UMI_DELIMITER)
        umi3_delimiter = ""
    else:
        umi5_ext_record = barcode_ext_record
    # Record after 3' UMI extraction.
    # Add UMI to record ID, using "_" delimiter for consistency
    # with UMI-tools, unless 5' UMI has been extracted, in which case
    # use "".
    umi3_ext_record = trim_fastq_record_3prime(
        umi5_ext_record,
        len(umi3),
        True,
        umi3_delimiter)
    return (record, trim_record, umi3_ext_record)


def create_fastq_examples(output_dir):
    """
    Create simulated fastq files.

    :param output_dir: Output directory
    :type output_dir: str or unicode
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    make_fastq_record("SRR", "AAAA")
    seed(42)  # Fix random seed so can repeatedly create same files

    # Components for simulated reads compatible with the vignette
    # files yeast_yAL_CDS_w_250utrs.fa.
    # These are aimed at the Duncan & Mata format with 4nt UMI at
    # each end of read.
    read_a = "ATGGCATCCACCGATTTCTCCAAGATTGAA"  # 30nt starting ORF of YAL003W
    read_ae = "ATGGCATCCACCGATGTCTCCAAGATTGAA"  # 1 error in read A
    read_b = "TCTAGATTAGAAAGATTGACCTCATTAA"  # 28nt immediately following start of ORF of YAL038W

    # UMIs
    umi_x = "CGTA"
    umi_y = "ATAT"
    umi_ye = "ATAA"  # 1 error in umi_y
    umi_z = "CGGC"
    umi_ze = "CTGC"  # 1 error in umi_x
    umi_5 = "AAAA"
    umi_5c = "CCCC"

    adaptor = "CTGTAGGCACC"  # Adaptor sequence used in vignette data

    post_adaptor_nt = "AC"

    # Create raw data with 5' and 3' UMIs and an adaptor.
    # M.N in the record names note the group the record is expected
    # belong to (M) and its number within that group.
    # After deduplication there should only be 1 member of each
    # group.
    config_5_3_adaptor = [
        ["EWSim-1.1-umi5-reada-umix", umi_5, read_a, umi_x, QUALITY_HIGH],
        ["EWSim-1.2-umi5-reada-umix", umi_5, read_a, umi_x, QUALITY_MEDIUM],
        ["EWSim-1.3-umi5-readae-umix", umi_5, read_ae, umi_x, QUALITY_MEDIUM],
        ["EWSim-2.1-umi5-reada-umiy", umi_5, read_a, umi_y, QUALITY_MEDIUM]
    ]
    # Create raw data with 5' and 3' UMIs and an adaptor plus an
    # extra nt past the "", adaptor for the shorter read.
    config_5_3_post_adaptor_nt = [
        ["EWSim-3.1-umi5-readb-umix", umi_5, read_b, umi_x, QUALITY_MEDIUM],
        ["EWSim-4.1-umi5-readb-umiz", umi_5, read_b, umi_z, QUALITY_HIGH],
        ["EWSim-4.2-umi5-readb-umiz", umi_5, read_b, umi_z, QUALITY_MEDIUM],
        ["EWSim-4.3-umi5-readb-umize", umi_5, read_b, umi_ze, QUALITY_MEDIUM],
        ["EWSim-5.1-umi5c-readb-umix", umi_5c, read_b, umi_x, QUALITY_MEDIUM]
    ]
    records = [
        make_fastq_records(tag, read, qualities, umi5, umi3, "", adaptor)
        for [tag, umi5, read, umi3, qualities] in config_5_3_adaptor]
    records_post_adaptor_nt = [
        make_fastq_records(tag, read, qualities, umi5, umi3, "",
                           adaptor, post_adaptor_nt)
        for [tag, umi5, read, umi3, qualities] in config_5_3_post_adaptor_nt]
    records.extend(records_post_adaptor_nt)
    file_names = ["example_umi5_umi3_umi_adaptor.fastq",
                  "example_umi5_umi3_umi.fastq",
                  "example_umi5_umi3.fastq"]
    for file_name, fastq_records in zip(file_names, zip(*records)):
        with open(os.path.join(output_dir, file_name), "w") as f:
            SeqIO.write(fastq_records, f, FASTQ_FORMAT)

    # Simulate raw data with only 3' umi.
    config_3 = [
        ["EWSim-1.1-reada-umix", read_a, umi_x, QUALITY_HIGH],
        ["EWSim-1.2-reada-umix", read_a, umi_x, QUALITY_MEDIUM],
        ["EWSim-1.3-readae-umix", read_ae, umi_x, QUALITY_MEDIUM],
        ["EWSim-2.1-reada-umiy", read_a, umi_y, QUALITY_MEDIUM],
        ["EWSim-3.1-readb-umix", read_b, umi_x, QUALITY_MEDIUM],
        ["EWSim-4.1-readb-umiz", read_b, umi_z, QUALITY_HIGH],
        ["EWSim-4.2-readb-umiz", read_b, umi_z, QUALITY_MEDIUM],
        ["EWSim-4.3-readb-umize", read_b, umi_ze, QUALITY_MEDIUM],
    ]
    records = [
        make_fastq_records(tag, read, qualities, "", umi3, "", adaptor)
        for [tag, read, umi3, qualities] in config_3]
    file_names = ["example_umi3_umi_adaptor.fastq",
                  "example_umi3_umi.fastq",
                  "example_umi3.fastq"]
    for file_name, fastq_records in zip(file_names, zip(*records)):
        with open(os.path.join(output_dir, file_name), "w") as f:
            SeqIO.write(fastq_records, f, FASTQ_FORMAT)

    # Create multiplexed data.
    # Use same data as 5' and 3' UMIs and an adaptor but with
    # barcodes.
    # Barcodes (keys) each with list of barcodes with 1-nt and 2-nt
    # mismatches.
    barcodes = {'ACG': ['ACT', 'TAG'],
                'GAC': ['GTC', 'GTA'],
                'CGA': ['TGA', 'CTT']}
    barcode_names = list(barcodes.keys())
    # Barcode that will be unassigned during demultiplexing.
    barcodes['TTT'] = []

    barcode_format = "-bar{:01d}.{:01d}"

    for file_name in ["example_multiplex_umi_barcode_adaptor.fastq",
                      "example_multiplex_umi_barcode.fastq",
                      "example_multiplex.fastq",
                      "example_multiplex_barcodes.tsv"]:
        file_path = os.path.join(output_dir, file_name)
        if os.path.exists(file_path):
            os.remove(file_path)
    with open(os.path.join(output_dir,
                           "example_multiplex_barcodes.tsv"), "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["SampleID", "TagRead"])
        for index, barcode in enumerate(barcode_names):
            writer.writerow(["Tag{:01d}".format(index), barcode])

    # Iterate over mismatches then barcodes so can interleave reads
    # for each barcode i.e. reads for each barcodes will be created
    # first then the reads for the 1nt mismatches then those for 2nt
    # mismatches.
    barcode_sets = [['ACG', 'GAC', 'CGA'],  # Barcodes
                    ['ACT', 'GTC', 'TGA'],  # 1nt mismatches
                    ['TAG', 'GTA', 'CTT']]  # 2nt mismatches
    for mismatch_index, barcode in enumerate(barcode_sets):
        for barcode_index, barcode in enumerate(barcodes):
            records = [
                make_fastq_records(tag +
                                   barcode_format.format(barcode_index,
                                                         mismatch_index),
                                   read, qualities,
                                   umi5, umi3, barcode,
                                   adaptor, "")
                for [tag, umi5, read, umi3, qualities] in config_5_3_adaptor]
            records_post_adaptor_nt = [
                make_fastq_records(tag +
                                   barcode_format.format(barcode_index,
                                                         mismatch_index),
                                   read, qualities,
                                   umi5, umi3, barcode,
                                   adaptor, post_adaptor_nt)
                for [tag, umi5, read, umi3, qualities] in config_5_3_post_adaptor_nt]
            records.extend(records_post_adaptor_nt)
            file_names = ["example_multiplex_umi_barcode_adaptor.fastq",
                          "example_multiplex_umi_barcode.fastq",
                          "example_multiplex.fastq"]
            for file_name, fastq_records in zip(file_names, zip(*records)):
                with open(os.path.join(output_dir, file_name), "a") as f:
                    SeqIO.write(fastq_records, f, FASTQ_FORMAT)

    # TODO GZIP


if __name__ == "__main__":
    create_fastq_examples(sys.argv[1])
