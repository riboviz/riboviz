# Developing and running tests

## Run vignette regression test suite

`riboviz.test.regression.test_vignette` is a regression test suite. The regression test suite runs `riboviz.tools.prep_riboviz` on the vignette data, in `vignette/`, then compares the results, in `vignette/`, to a directory of pre-calculated results, specified by the user.

Usage:

```console
$ pytest riboviz/test/regression/test_vignette.py \
    --expected=DIRECTORY \
    [--skip-workflow] \
    [--check-index-tmp]
```

The test suite accepts three custom command-line parameters:

* `--expected=<DIRECTORY>`: Directory with expected data files, against which files in vignette/ will be checked.
* `--skip-workflow`: Workflow will not be run prior to checking data files.
* `--check-index-tmp`: Check index and temporary files (default is that only the output files are checked).

If `--skip-workflow` is not specified then `riboviz.tools.prep_riboviz` is run using the vignette configuration, `vignette/vignette-config.yaml`.

The vignette output files (and the index and temporary files, if `--check-index-tmp` was provided) in `vignette/` are then compared against those in the directory provided via the expected parameter.

If `--check-index-tmp` is not provided (the default behaviour) then tests for index and temporary files will be skipped. This will appear as follows:

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

### Regression test data

Regression test data can be found within repositories prefixed `regression-test-data-<YYYYMMDD>`, within the [riboviz](https://github.com/riboviz) project on GitHub. These contain output files from the workflow.

A regression test data repository can be used as follows (for example):

```console
$ git clone https://github.com/riboviz/regression-test-data-<YYYYMMDD>
$ cd riboviz
$ pytest riboviz/test/regression/test_vignette.py --expected=$HOME/regression-test-data-<YYYYMMDD>
```

**Note:** The most recent `regression-test-data-<YYYYMMDD>` repository will be consistent with the most recent version of the code on the `develop` branch.

**Note:** Do not use the `--check-index-tmp` parameter when running regression tests using these repositories as the repositories contain no index or temporary files.

### Using your own regression test directory

After running the vignette you can copy `vignette` directory and then use that directory as the argument provided to the `--expected` parameter. This can be useful if you are changing any aspect of the workflow and want to ensure that you have not inadvertently changed the contents of the index or temporary files in any way.

For example:

```console
$ cp -r vignette/ ~/regression-test-data-<BRANCH>-<COMMIT-HASH>
...after doing some development...
$ pytest riboviz/test/regression/test_vignette.py \
    --expected=$HOME/regression-test-data-<BRANCH>-<COMMIT-HASH> \
    --check-index-tmp
```

Alternatively, you may find this useful if your changes *were* intended to change the contents of the index, temporary, or output, files - the copy of `vignette` serves as your new regression test directory. You could use this as a basis for [Creating a regression test data repository](./create-test-data-repository.md)

---

## Run all tests (excluding regression tests)

Run:

```console
$ pytest -v --ignore-glob="*regression*"
```

Run all tests and create a test coverage report which includes the line numbers of statements that were not executed:

```console
$ pytest --cov-config=.coveragerc --cov-report term-missing \
    --cov=riboviz --ignore-glob="*regression*"
```

---

## Useful pytest flags

* `-s`: disable output capture so, for example, `print` messages are shown.
* `-v`: verbose mode, displays names of test functions run.
* `-k`: run a specific test function.
