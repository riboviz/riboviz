# Developing and running tests

## Run vignette regression test suite

`riboviz.test.regression.test_vignette` is a regression test suite. The regression test suite runs riboviz.tools.prep_riboviz on the vignette data, in `vignette/`, then compares the results, in `vignette/`, to a directory of pre-calculated results, specified by the user.

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
