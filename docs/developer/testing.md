# Developing and running tests


* [Run Python tests and workflow tests](#run-python-tests-and-workflow-tests)
* [Run R tests](#run-r-tests)
* [Download vignette integration test data](#download-vignette-integration-test-data)
* [Run vignette integration tests](#run-vignette-integration-tests)
* [Using the integration test suite](#using-the-integration-test-suite)
  - [Specifying values for environment variable configuration tokens](#specifying-values-for-environment-variable-configuration-tokens)
- [Skipping tests for index and temporary files](#skipping-tests-for-index-and-temporary-files)
  - [Using your own expected results directory](#using-your-own-expected-results-directory)
  - [How actual directories and files are compared to expected directories and files](#how-actual-directories-and-files-are-compared-to-expected-directories-and-files)
  - [Limitations of tests for UMI extraction, deduplication and grouping](#limitations-of-tests-for-umi-extraction-deduplication-and-grouping)
* [Useful pytest flags](#useful-pytest-flags)

---

## Run Python tests and workflow tests

Run:

```console
$ pytest --ignore-glob="*integration*"
```

`PendingDeprecationWarning` `warnings` can be ignored.

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

---

## Download vignette integration test data

Integration test data can be found within repositories prefixed `test-data-<YYYYMMDD>|<TAG>`, within the [riboviz](https://github.com/riboviz) project on GitHub. These contain output files from the workflow.

**Note:** The most recent `test-data-<YYYYMMDD>` repository will be consistent with the most recent version of the code on the `develop` branch.

**Note:** `test-data-<TAG>` repositories will be consistent with releases with tagged with `<TAG>`.

**Note:** Repositories prior to release 2.1, June 2021, were prefixed `regression-`.

Download test data, for example:

```console
$ cd
$ git clone https://github.com/riboviz/test-data-2.1
```

---

## Run vignette integration tests

These assume you have vignette integration test data. See [Download vignette integration test data](#download-vignette-integration-test-data) above.

Run the integration tests for the workflow (these may take a few minutes):

```console
$ cd riboviz
$ pytest riboviz/test/integration/test_integration.py \
    --expected=$HOME/test-data-<YYYYMMDD>|<TAG> \
    --config-file=vignette/vignette_config.yaml
```

For example:

```console
$ cd riboviz
$ pytest riboviz/test/integration/test_integration.py \
    --expected=$HOME/test-data-2.1 \
    --config-file=vignette/vignette_config.yaml
```

**Note:** If you have already run the vignette, then you can add a `--skip-workflow` flag to the call to `pytest`.

**Troubleshooting: `PendingDeprecationWarning`**

`PendingDeprecationWarning` `warnings` can be ignored.

**Troubleshooting: `FAILED ... test_bam_to_h5_h5[...]` test failures**

If, running the vignette integration tests, you see the following two failures:

```
...
FAILED riboviz/test/integration/test_integration.py::test_bam_to_h5_h5[vignette/output-WTnone]
FAILED riboviz/test/integration/test_integration.py::test_bam_to_h5_h5[vignette/output-WT3AT]
...
```

but all other tests pass, then these two test failures can be ignored. These failures can arise if the test data was produced in an environment using rhdf 2.34.0 or above and you are using a version of rhd5 prior to 2.34.0, or vice versa. [Bioconductor 3.12 Released](http://bioconductor.org/news/bioc_3_12_release/) (October 28, 2020) explains that, for rhdf4 2.34.0, "datasets written with h5write() now have the attribute rhdf5-NA.OK added to them ... to indicate that rhdf5 was used to create the file...". These new attributes have form:

```
ATTRIBUTE "rhdf5-NA.OK" {
   DATATYPE  H5T_STD_I32LE
   DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
   DATA {
   (0): 1
   }
}
```

The H5 files produced by riboviz (via invocation of the `bam_to_h5.R` script) are used by downstream processing steps (notably `generate_stats_figs.R`) so these failures can be ignored if all other integration tests pass.

---

## Using the integration test suite

`riboviz.test.integration.test_integration` is a integration test suite. The test suite runs the workflow using a given configuration file, then compares the results to a directory of pre-calculated results, specified by the user.

Usage:

```console
$ pytest riboviz/test/integration/test_integration.py \
    --expected=<EXPECTED_RESULTS_DIRECTORY> \
    [--skip-workflow] \
    [--check-index-tmp] \
    [--config-file=<CONFIG_FILE>]
```

The test suite accepts the following command-line parameters:

* `--expected=<EXPECTED_RESULTS_DIRECTORY>`: Directory with expected data files, against which files specified in the configuration file (see below) will be checked.
* `--skip-workflow`: Workflow will not be run prior to checking data files. This can be used to check existing files generated by a run of the workflow.
* `--check-index-tmp`: Check index and temporary files (default is that only the output files are checked).
* `--config-file=<CONFIG_FILE>`: Configuration file. If provided then the index, temporary and output directories specified in this file will be validated against those specified by `--expected`. If not provided then the file `vignette/vignette_config.yaml` will be used.

**Note:** Do not use the `--check-index-tmp` parameter when running the tests using the test data repositories `https://github.com/riboviz/test-data-<YYYYMMDD>|<TAG>`. These repositories contain no index or temporary files.

### Specifying values for environment variable configuration tokens

To specify values for environment variables cited as tokens in configuration parameters (see [Environment variables and configuration tokens](../user/prep-riboviz-config.md#environment-variables-and-configuration-tokens)), then these should be defined as described in [Defining values for environment variables](../user/prep-riboviz-run-nextflow.md#defining-values-for-environment-variables). For example:

```console
$ RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY> \
  RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY> \
  RIBOVIZ_DATA=<DATA_DIRECTORY> \
  pytest riboviz/test/integration/test_integration.py \
    --expected=<EXPECTED_RESULTS_DIRECTORY> \
    [--skip-workflow] \
    [--check-index-tmp] \
    [--config-file=<CONFIG_FILE>]
```

### Skipping tests for index and temporary files

If `--check-index-tmp` is not provided (the default behaviour) then tests for index and temporary files will be skipped. An example of how this appears is as follows:

```
...
riboviz/test/integration/test_integration.py::test_hisat2_build_index[vignette/index-True-YAL_CDS_w_250-1] SKIPPED
...
riboviz/test/integration/test_integration.py::test_hisat2_build_index[vignette/index-True-yeast_rRNA-8] SKIPPED
riboviz/test/integration/test_integration.py::test_cutadapt_fq[vignette/tmp-False-WTnone] SKIPPED
...
riboviz/test/integration/test_integration.py::test_samtools_index_dedup_bam[vignette/tmp-False-WT3AT] SKIPPED
riboviz/test/integration/test_integration.py::test_samtools_view_sort_index[vignette/output-False-WTnone] PASSED
riboviz/test/integration/test_integration.py::test_samtools_view_sort_index[vignette/output-False-WT3AT] PASSED
riboviz/test/integration/test_integration.py::test_umitools_dedup_stats_tsv[vignette/tmp-False-False-WTnone-edit_distance.tsv] SKIPPED
...
```

### Using your own expected results directory

After running a workflow using your own configuration file you can copy your index, temporary and output directories and then use those copies for future tests. One way to do this is to create a new folder with these results e.g.

```console
$ mkdir results-<BRANCH>-<COMMIT-HASH>
$ cp -r <INDEX_DIRECTORY> results-<BRANCH>-<COMMIT-HASH>/
$ cp -r <TMP_DIRECTORY> results-<BRANCH>-<COMMIT-HASH>/
$ cp -r <OUTPUT_DIRECTORY> results-<BRANCH>-<COMMIT-HASH>/
```

**Note:** Keep the names of the original index, temporary and output directories. See [How actual directories and files are compared to expected directories and files](#how-actual-directories-and-files-are-compared-to-expected-directories-and-files) above).

You can now provide this folder as a value for the `--expected` parameter. For example, after doing some development, you can run:

```console
$ pytest riboviz/test/integration/test_integration.py \
    --expected=results-<BRANCH>-<COMMIT-HASH> \
    --config-file=<CONFIG_FILE> \
    --check-index-tmp
```

### How actual directories and files are compared to expected directories and files

As the expected data directories (`--expected`) and those with the data to be tested may vary in their paths the following approach is used:

* The paths of the directories with the data to be tested are taken to be those specified in the configuration file.
* The paths of the directories with the expected data are taken to be relative to the `--expected` directory and to share common names with the final directory names in each path of the actual data directories.

For example, if the configuration file has:

```yaml
dir_index: vignette/index
dir_out: vignette/output
dir_tmp: vignette/tmp
```

and `--expected` is `$HOME/results-mybranch-abcdefg` then the directories with the data to be tested are:

```
vignette/index
vignette/output
vignette/tmp
```

and the directories with the expected data are:

```
$HOME/results-mybranch-abcdefg/index
$HOME/results-mybranch-abcdefg/output
$HOME/results-mybranch-abcdefg/tmp
```

Observe that the final directories in each path - `index`, `output`, `tmp` - are the same, though the paths to these differ.

### Limitations of tests for UMI extraction, deduplication and grouping

If running with a configuration that used UMI extraction, deduplication and grouping then:

* UMI deduplication statistics files are not checked (files prefixed by `dedup_stats`).
* UMI group file post-deduplication files, (`post_dedup_groups.tsv`) files can differ between runs depending on which reads are removed by `umi_tools dedup`, so only the existence of the file is checked.
* BAM file output by deduplication (`dedup.bam`) files can differ between runs depending on which reads are removed by `umi_tools dedup`, so only the existence of the file is checked.

---

## Useful pytest flags

* `-s`: disable output capture so, for example, `print` messages are shown.
* `-v`: verbose mode, displays names of test functions run.
* `-k`: run a specific test function.
* `--cov-config=.coveragerc --cov-report term-missing --cov=riboviz`: create a test coverage report which includes the line numbers of statements that were not executed.
