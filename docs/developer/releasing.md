# Building a release

These are the key steps in building a release.

* Merge all pull requests agreed for the release into the `develop` branch, via the GitHub web interface.
* Merge additional branches agreed for the release into the `develop` branch via creation of pull requests, again via the GitHub web interface.
* Check [docs/acks.md](../acks.md), and update if necessary.
* Check [docs/reference/references.md](../reference/references.md), and update if necessary.
* Update the following documentation to refer to the current release version:
  - [docs/developer/create-test-data-repository.md](./create-test-data-repository.md)
  - [docs/developer/testing.md](./testing.md)
  - [docs/developer/releasing.md](./releasing.md) (this page)
* Update the Releases table in `README.md` with a row for intended release and provide a URL with the intended release tag. See the current rows of the table for the format of the URL.
* Update software versions in `docs/user/install.md` - see [Update Dependencies overview in Install RiboViz and dependencies](./documentation.md#update-dependencies-overview-in-install-riboviz-and-dependencies)
* Merge `develop` branch into `master` branch.
* Rerun all tests:
  - [Run Python tests and workflow tests](./testing.md#run-python-tests-and-workflow-tests).
  - [Run R tests](./testing.md#run-r-tests).
  - [Run vignette regression tests](./testing.md#run-vignette-regression-tests).
* Skim all documentation via the GitHub web interface and check formatting and links.
* Tag release in repository:

```console
$ git tag 2.0
$ git push origin 2.0
```

* Create a regression test data repository. See [Creating a regression test data repository](./create-test-data-repository.md)
