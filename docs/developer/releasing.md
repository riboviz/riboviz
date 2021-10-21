# Building a release

The process for building a release is as follows:

1. [Merge agreed pull requests and branches](#1-merge-agreed-pull-requests-and-branches)
2. [Check and update documentation](#2-check-and-update-documentation)
3. [Test](#3-test)
4. [Release](#4-release)

---

## 1. Merge agreed pull requests and branches

* Merge all pull requests agreed for the release into the `develop` branch, via the GitHub web interface.
* Merge additional branches agreed for the release into the `develop` branch via creation of pull requests, again via the GitHub web interface.

---

## 2. Check and update documentation

* Update `README.md`:
  - Update the 'figshare. Software' reference publication year and authors list, if required.
  - Update the 'Copyright' year, if required.
  - Update the 'Releases' table with a row for the release and provide a URL with the release tag. See the current table rows for the URL format.
* Update documentation to refer to the current release version:
  - `docs/developer/releasing.md` (this page)
  - `docs/developer/create-test-data-repository.md`
  - `docs/developer/integration-tests.md`. Update the test data URLs in the `git clone` commands with the intended URL for the test data repository (`https://github.com/riboviz/test-data-<TAG>`).
* Update acknowledgements and references, where necessary:
  - `docs/acks.md`
  - `docs/reference/references.md`
* Update software versions in `docs/user/install.md`. Ssee [Updating Dependencies overview tables](./dependencies.md#updating-dependencies-overview-tables) in [Adding and updating dependencies](./dependencies.md).
* Skim all documentation via the GitHub web interface and check, and fix, formatting and links.

---

## 3. Test

* [Run Python tests and workflow tests](./dev-python.md#run-python-tests-and-workflow-tests).
* [Run R tests](./dev-r.md#run-r-tests).
* [Run vignette integration tests](./integration-tests.md#run-vignette-integration-tests).

---

## 4. Release

* Create a test data repository. See [Creating a test data repository](./create-test-data-repository.md). Ensure this has the same name as that in [docs/developer/create-test-data-repository.md](./create-test-data-repository.md) which you updated above.
* Merge `develop` branch into `main` branch.
* Tag release in repository:

```console
$ git tag 2.1
$ git push origin 2.1
```

* Create release on GitHub:
  - Visit https://github.com/riboviz/riboviz
  - Click Releases, https://github.com/riboviz/riboviz/releases
  - Click Draft a new release
  - Enter Tag version: 2.1
  - Enter Release title: 2.1
  - Enter Description: riboviz release 2.1 (RELEASE-DATE).
  - Click Publish release.
