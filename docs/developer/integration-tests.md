# Developing and running integration tests

* [Download vignette integration test data](#download-vignette-integration-test-data)
* [Run vignette integration tests](#run-vignette-integration-tests)
* [Using the integration test suite](#using-the-integration-test-suite)
  - [Specifying values for environment variable configuration tokens](#specifying-values-for-environment-variable-configuration-tokens)
  - [Skipping tests for index and temporary files](#skipping-tests-for-index-and-temporary-files)
  - [Using your own expected results directory](#using-your-own-expected-results-directory)
  - [How actual directories and files are compared to expected directories and files](#how-actual-directories-and-files-are-compared-to-expected-directories-and-files)
  - [Limitations of tests for UMI extraction, deduplication and grouping](#limitations-of-tests-for-umi-extraction-deduplication-and-grouping)
* [Writing an integration test](#writing-an-integration-test)
  - [Integration test fixtures](#integration-test-fixtures)
  - [Integration test parameters](#integration-test-parameters)
  - [Comparing expected and actual files for equality](#comparing-expected-and-actual-files-for-equality)
  - [Anatomy of an integration test function](#anatomy-of-an-integration-test-function)
  - [Conditionally skipping tests](#conditionally-skipping-tests)
* [Useful pytest flags](#useful-pytest-flags)
* [Information on pytest fixtures and parameters](#information-on-pytest-fixtures-and-parameters)

---

## Download vignette integration test data

Integration test data can be found within repositories prefixed `test-data-<YYYYMMDD>|<TAG>`, within the [riboviz](https://github.com/riboviz) project on GitHub. These contain output files from the workflow.

**Note:** The most recent `test-data-<YYYYMMDD>` repository will be consistent with the most recent version of the code on the `develop` branch.

**Note:** `test-data-<TAG>` repositories will be consistent with releases with tagged with `<TAG>`.

**Note:** Repositories prior to release 2.1, June 2021, were prefixed `regression-`.

Download test data, for example:

```console
$ cd
$ git clone https://github.com/riboviz/test-data-2.2
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
    --expected=$HOME/test-data-2.2 \
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

`riboviz.test.integration.test_integration` is an integration test suite. The test suite runs the workflow using a given configuration file, then compares the results to a directory of pre-calculated results, specified by the user.

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

After running a workflow using your own configuration file you can copy your index, temporary and output directories and then use those copies for future tests. One way to do this is to create a new folder with these results. For example:

```console
$ mkdir results-<BRANCH>-<COMMIT-HASH>
$ cp -r <INDEX_DIRECTORY> results-<BRANCH>-<COMMIT-HASH>/
$ cp -r <TMP_DIRECTORY> results-<BRANCH>-<COMMIT-HASH>/
$ cp -r <OUTPUT_DIRECTORY> results-<BRANCH>-<COMMIT-HASH>/
```

For example:

```console
$ mkdir results-develop-abcdefg
$ cp -r vignette/index results-develop-abcdefg/
$ cp -r vignette/tmp results-develop-abcdefg/
$ cp -r vignette/output results-develop-abcdefg/
```

**Note:** Keep the names of the original index, temporary and output directories. See [How actual directories and files are compared to expected directories and files](#how-actual-directories-and-files-are-compared-to-expected-directories-and-files) below).

You can now provide this folder as a value for the `--expected` parameter. For example, after doing some development, you can run:

```console
$ pytest riboviz/test/integration/test_integration.py \
    --expected=results-<BRANCH>-<COMMIT-HASH> \
    --config-file=<CONFIG_FILE> \
    --check-index-tmp
```

For example:

```console
$ pytest riboviz/test/integration/test_integration.py \
    --expected=results-develop-abcdefg \
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

and `--expected` is `$HOME/results-develop-abcdefg` then the directories with the data to be tested are:

```
vignette/index
vignette/output
vignette/tmp
```

and the directories with the expected data are:

```
$HOME/results-develop-abcdefg/index
$HOME/results-develop-abcdefg/output
$HOME/results-develop-abcdefg/tmp
```

Observe that the final directories in each path - `index`, `output`, `tmp` - are the same, though the paths to these differ.

### Limitations of tests for UMI extraction, deduplication and grouping

If running with a configuration that used UMI extraction, deduplication and grouping then:

* UMI deduplication statistics files are not checked (files prefixed by `dedup_stats`).
* UMI group file post-deduplication files, (`post_dedup_groups.tsv`) files can differ between runs depending on which reads are removed by `umi_tools dedup`, so only the existence of the file is checked.
* BAM file output by deduplication (`dedup.bam`) files can differ between runs depending on which reads are removed by `umi_tools dedup`, so only the existence of the file is checked.

---

## Writing an integration test

Here, the key parts of the integration tests are described. This, along with the code already in `riboviz/test/integration/test_integration.py`, is intended to provide the information required to edit existing or write new integration tests to validate temporary or output files by comparing these to temporary or output files in an integration test data directory.

### Integration test fixtures

All integration test functions must use the fixture, `prep_riboviz_fixture`, (defined in `riboviz/test/integration/test_integration.py`) which ensures the workflow is run if the `--skip-workflow` command-line parameter is not provided when running the integration tests. This is a module-wide fixture so is run once per invocation of `riboviz.test.integration.test_integration`. The fixture should be specified before the test function declaration as follows:

```python
@pytest.mark.usefixtures("prep_riboviz_fixture")
```

All integration test functions to validate temporary files must use the fixture, `skip_index_tmp_fixture` (defined in `riboviz/test/integration/conftest.py`) which ensures the test is skipped if the `--check-index-tmp` command-line parameter is provided when running the integration tests. This is a module-wide fixture so is run once per invocation of `riboviz.test.integration.test_integration`. The fixture should be specified before the test function declaration as follows:

```python
@pytest.mark.usefixtures("skip_index_tmp_fixture")
```

These fixtures are defined before the test functions as the test functions do not need to rely on values provided by the fixtures, only on their side-effects.

Additional fixtures, from which integration test functions can take values by declaring arguments with the same name as the fixture, are as follows:

| Fixture name | Description | Scope | Definition file |
| ------------ | ----------- | ----- | --------------- |
| `expected_fixture` | Value of `--expected` command-line option when the integration tests are run i.e., the integration test data. | module | `riboviz/test/integration/conftest.py` |
| `config_fixture` | Value of `--config-file` command-line option when the integration tests are run (default `vignette/vignette_config.yaml`). | module | `riboviz/test/integration/conftest.py` |
| `scratch_directory` | Scratch directory, created as a sub-directory of `tmpdir`, see below. | function | `riboviz/test/integration/test_integration.py` |
| `tmpdir` | Temporary directory, unique to test invocation. | function | Provided by `pytest`, see [Temporary directories and files](https://docs.pytest.org/en/6.2.x/tmpdir.html).

### Integration test parameters

The following parameters are available to parameterise integration test functions. These are derived from the YAML configuration file provided by the `--config-file` command-line option when the integration tests are run (default `vignette/vignette_config.yaml`). Test functions can use these parameters by declaring arguments with the same name as the parameter. The parameters (defined in `riboviz/test/integration/conftest.py`) are as follows.

| Fixture name | Description |
| ------------ | ----------- |
| `sample` | If `fq_files` is defined in the configuration, then this parameter parameter has the sample names from this value. Else if `multiplex_fq_files` is defined in the configuration then this parameter has the sample names that are deduced from the names of directories in `dir_out`, within the integration test data directory, cross-referenced with the sample sheet file specified in `sample_sheet`. If any sample name is `NotHere` then it is removed. A test taking this parameter will be executed once for each sample in turn. |
| `is_multiplexed` | `True` if `multiplex_fq_files` in the configuration defines one or more files, `False` otherwise. |
| `multiplex_name` | Multiplexed file names prefixes, without extensions, from `multiplex_fq_files`, if any. A test taking this parameter will be executed for each such file in turn. |
| `index_prefix` | Indexed file prefix values (`orf_index_prefix` and `rrna_index_prefix`). A test taking this parameter will be executed for each prefix in turn. |
| `<param>` | `<param>` is a configuration parameter name. Its value will be taken from the configuration. For undefined values, default values are taken from `riboviz/default_config.yaml`. |

### Comparing expected and actual files for equality

The following functions are available for comparing files. See the function definitions for more information on the nature of the comparisons done.

| Definition File | Function |
| --------------- | -------- |
| `riboviz/bedgraph.py` | `equal_bedgraph(file1, file2)` |
| `riboviz/count_reads.py` | `equal_read_counts(file1, file2, comment="#")` |
| `riboviz/fastq.py` | `equal_fastq(file1, file2)` |
| `riboviz/h5.py` | `equal_h5(file1, file2)` |
| `riboviz/html.py` | `equal_html(file1, file2)` |
| `riboviz/sam_bam.py` | `equal_bam(file1, file2)` |
| `riboviz/sam_bam.py` | `equal_sam(file1, file2)` |
| `riboviz/utils.py` | `equal_file_names(file1, file2)` |
| `riboviz/utils.py` | `equal_file_sizes(file1, file2)` |
| `riboviz/utils.py` | `equal_tsv(file1, file2, tolerance=0.0001, ignore_row_order=False, comment="#", na_to_empty_str=False)` |

Complementing these, the following helper functions are defined in `riboviz/test/integration/test_integration` for comparing files. See the function definitions for more information on the nature of the comparisons done.

```python
compare_tsv_files(expected_fixture, directory, subdirectory, file_name)
compare_fq_files(expected_fixture, directory, subdirectory, file_name)
compare_sam_files(expected_directory, directory, scratch_directory, sample, file_name)
check_pdf_file_exists(dir_out, sample, file_name)
```

### Anatomy of an integration test function

As an example consider the following test:

```python
@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name", [
    workflow_files.ORF_MAP_SAM,
    workflow_files.RRNA_MAP_SAM])
def test_hisat2_sam(expected_fixture, dir_tmp, scratch_directory,
                    sample, file_name):
    compare_sam_files(expected_fixture, dir_tmp, scratch_directory,
                      sample, file_name)
```

This can be broken down as follows.

```python
@pytest.mark.usefixtures("skip_index_tmp_fixture")
```

The test function validates temporary files, so it uses the fixture `skip_index_tmp_fixture`, so it can be skipped if the `--check-index-tmp` command-line parameter is provided when running the integration tests.

```python
@pytest.mark.usefixtures("prep_riboviz_fixture")
```

In common with all integration test functions, it uses the fixture `prep_riboviz_fixture`, so that the workflow will be run prior to running the tests, unless the `--skip-workflow` command-line parameter is provided when running the integration tests.

```python
@pytest.mark.parametrize("file_name", [
    workflow_files.ORF_MAP_SAM,
    workflow_files.RRNA_MAP_SAM])
def test_hisat2_sam(expected_fixture, dir_tmp, scratch_directory,
                    sample, file_name):
```

The test function is parameterised to take two file names i.e. it will run twice, the first time with `workflow_files_ORF_MAP_SAM` (which has value `orf_map.sam`), the second time with `workflow_files.RRNA_MAP_SAM` (which has value `rRNA_map.sam`). It also takes the following fixtures and parameters:

* `expected_fixture`: fixture providing the location of the integration test data against which files are to be validated.
* `dir_tmp`: parameter providing the value of the `dir_tmp` configuration parameter.
* `scratch_directory`: fixture providing a scratch directory, a sub-directory of a temporary directory created by pytest.
* `sample`: parameter providing the name of each sample in turn.

If there are three samples defined in `fq_files` e.g., `WTnone` and `WT3AT`, and as `file_name` has values `orf_map.sam` `rRNA_map.sam` then the combination of these parameters means that `test_hisat2_sam`  would be run for each of the following combinations of parameters

| `sample` | `file_name`    | 
| -------- | -------------- |
| `WTnone` | `orf_map_sam`  |
| `WTnone` | `rRNA_map_sam` |
| `WT3AT`  | `orf_map_sam`  |
| `WT3AT`  | `rRNA_map_sam` |

### Conditionally skipping tests

Tests can be conditionally skipped, based on the value of configuration parameters, as follows:

* Declare the configuration parameter as an argument to the test function. It will then be passed in as a parameter, as described in [Integration test parameters](#integration-test-parameters) above.
* Declare a conditional invocation of `pytest.skip` to tell pytest to skip the test depending on the desired condition.

For example, the file `normalized_density_APEsites_per_codon.pdf` is only output by the workflow (specifically `generate_stats_figs.R`) if values for the configuration parameters `t_rna_file` and `codon_positions_file` are provided and also if the configuration parameter `output_pdfs` is `True`. If any of these conditions do not hold then no output file is produced and so, correspondingly, the test for this file should be skipped. The implementation of the corresponding test is as follows:

```python
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize(
    "file_name",
    [workflow_r.NORMALIZED_DENSITY_APESITES_PER_CODON_PDF])
def test_generate_stats_figs_t_rna_codon_positions_pdf(
        t_rna_file, codon_positions_file, output_pdfs, dir_out,
        sample, file_name):
    if not t_rna_file:
        pytest.skip('Skipped test as t_rna_file: {}'.format(
            t_rna_file))
    if not codon_positions_file:
        pytest.skip('Skipped test as codon_positions_file: {}'.format(
            codon_positions_file))
    if not output_pdfs:
        pytest.skip('Skipped test as output_pdfs: {}'.format(output_pdfs))
    check_pdf_file_exists(dir_out, sample, file_name)
```

## Useful pytest flags

See [Useful pytest flags](./dev-python.md#useful-pytest-flags) in [Developing Python components](./dev-python.md).

---

## Information on pytest fixtures and parameters

See [Information on pytest fixtures and parameters](./dev-python.md#information-on-pytest-fixtures-and-parameters) in [Developing Python components](./dev-python.md).
