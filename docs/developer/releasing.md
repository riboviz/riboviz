# Building a release

These are the key steps in building a release.

* Merge all pull requests agreed for the release into the `develop` branch, via the GitHub web interface.
* Merge additional branches agreed for the release into the `develop` branch via creation of pull requests, again via the GitHub web interface.
* Check [docs/acks.md](../acks.md), and update if necessary.
* Check [docs/reference/references.md](../reference/references.md), and update if necessary.
* Update software versions in `docs/user/install.md` - see [Update Dependencies overview](./dependencies.md#update-dependencies-overview).
* Update the following documentation to refer to the current release version:
  - [docs/developer/releasing.md](./releasing.md) (this page)
  - [docs/developer/create-test-data-repository.md](./create-test-data-repository.md)
  - [docs/developer/testing.md](./testing.md), update the test data URLs in the `git clone` commands with the intended URL for the test data repository (`https://github.com/riboviz/test-data-<TAG>`).
* Update `README.md`:
  - Update the 'figshare. Software' reference publication year.
  - Update the 'Copyright' year.
  - Update the Releases table with a row for the release and provide a URL with the release tag. See the current rows of the table for the format of the URL.
* Skim all documentation via the GitHub web interface and check formatting and links.
* Rerun all tests:
  - [Run Python tests and workflow tests](./testing.md#run-python-tests-and-workflow-tests).
  - [Run R tests](./testing.md#run-r-tests).
  - [Run vignette integration tests](./testing.md#run-vignette-integration-tests).
* Create a test data repository. See [Creating a test data repository](./create-test-data-repository.md). Ensure this has the same name as that in [docs/developer/create-test-data-repository.md](./create-test-data-repository.md) which you updated above.
* Merge `develop` branch into `main` branch.
* Tag release in repository:

```console
$ git tag 2.1
$ git push origin 2.1
```
