# Building a release

These are the key steps in building a release.

* Merge all pull requests agreed for the release into the `develop` branch, via the GitHub web interface.
* Merge additional branches agreed for the release into the `develop` branch via creation of pull requests, again via the GitHub web interface.
* Merge `develop` branch into `master` branch.
* Update the Releases table in `README.md` with a row for intended release and provide a URL with the intended release tag. See the current rows of the table for the format of the URL.
* Update software versions in `docs/user/install.md` - see [Update Dependencies overview in Install RiboViz and dependencies](./documentation.md#update-dependencies-overview-in-install-riboviz-and-dependencies)
* Rerun all tests - see [Run all tests (excluding regression tests)](./testing.md#run-all-tests-excluding-regression-tests).
* Rerun Python workflow and regression tests - see [Run vignette regression test suite](./testing.md#run-vignette-regression-test-suite).
* Rerun Nextflow workflow and regression tests - see [Run vignette regression test suite](./testing.md#run-vignette-regression-test-suite).
* Skim all documentation via the GitHub web interface and check formatting and links.
* Tag release in repository:

```console
$ git tag 2.0
$ git push origin 2.0
```

* Create a regression test data repository. See [Creating a regression test data repository](./create-test-data-repository.md)
