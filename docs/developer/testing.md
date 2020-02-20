# Developing and running tests

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

## Run non-regression tests

Run all tests (excluding regression tests):

```console
$ pytest -v --ignore-glob="*regression*"
```

Run all tests and generate test coverage report (including lines that were not invoked):

```console
$ pytest --cov-config=.coveragerc --cov-report term-missing \
    --cov=riboviz --ignore-glob="*regression*"
```
