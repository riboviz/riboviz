# Building a release

The process for building a release is as follows:

1. [Merge agreed pull requests and branches](#1-merge-agreed-pull-requests-and-branches)
2. [Check and update documentation](#2-check-and-update-documentation)
3. [Test](#3-test)
4. [Release](#4-release)
5. [Update figshare deposit](#5-update-figshare-deposit)

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

---

## 5. Update figshare deposit

**Note:** This step can only be done by members of the private figshare 'riboviz: analysis and visualization of ribosome profiling datasets' project.

Download release from GitHub:

* Visit https://github.com/riboviz/riboviz/tags
* Click 'zip' for version to deposit.
* This downloads URL https://github.com/riboviz/riboviz/archive/refs/tags/TAG.zip into a file `riboviz-<TAG>.zip`.

Update riboviz article:

* Log in to [figshare](https://figshare.com/).
* Click Projects.
* Click 'riboviz: analysis and visualization of ribosome profiling datasets'.
* In the content table, 'riboviz: software for analysis and visualization of ribosome profiling datasets' row, click the cog icon and select Edit item.
* Click the 'X' icon to next to `riboviz-<OLD_TAG>.zip` to remove that file. Previous deposits will not be affected.
* Click browse.
* Select `riboviz-<TAG>.zip` that you downloaded above.
* Edit any fields that need updating e.g. Authors, Description, Funding, References etc.
* Double-check that titles, authors (including spellings) and files are correct before publication. Correcting these after publication results in a new version being created by figshare. For more information, see figshare's [Can I edit or delete my research after it has been made public?](https://help.figshare.com/article/can-i-edit-or-delete-my-research-after-it-has-been-made-public).
* Click Save changes.
* Check Publish changes.

**Note:** Each specific fighsare deposit has its own DOI e.g., [10.6084/m9.figshare.12624200.v1](https://doi.org/10.6084/m9.figshare.12624200.v1), [10.6084/m9.figshare.12624200.v2](https://doi.org/10.6084/m9.figshare.12624200.v2) etc, but the base fighare DOI, [10.6084/m9.figshare.12624200](https://doi.org/10.6084/m9.figshare.12624200), as listed in `README.md` remains unchanged and always points to the most recent deposit.
