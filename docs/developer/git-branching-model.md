# Git branching model

RiboViz uses three classes of branch for development:

* `master` branch: the most stable branch. Releases are created as tags of specific versions from `master`.
* `develop` branch: for ongoing development. At regular intervals this will be merged into `master`.
* Feature and bug fix-specific branches. On completion, these will be merged into `develop`.

## Developing new features, enhancements and bug fixes

When developing new features, enhancements and bug fixes, developers should create a new development branch from `develop`.

For branches associated with GitHub issues, the naming scheme `<name>-<#>` is used, where:

* `<name>` is a short name for the issue. This should be lower-case, with `-` used as a delimiter if desired.
* `<#>` is the number of the issue as in the GitHub issue tracker.
* If multiple developers are working on the same issue branch then a developer may want their own branch off this issue branch. These can be named `<name>-<#>-<user-name>`.

For example, for the issue "Investigate cutadapt -j flag #43" the branch name was `configure-cutadapt-cores-43`. If `mikej888` wanted his own branch off this issue branch, then it could be called `configure-cutadapt-cores-43-mikej888`.

Preferably, all branches should have an associated issue.

Please do **not** commit large (many megabytes or gigabytes) data files into the repository. Data files that are a few megabytes in size *are* permitted if required for test code.

Please be reassured that so long as you do **not** merge branches into `develop` or `master` or commit large data files then any mistakes you make should be relatively straightforward to recover from!

## Merging development branches into `develop`

To request that a branch be merged into `develop`:

1. Merge the current version of `develop` into your local branch.
2. Resolve conflicts, if any, and commit the fixes.
3. [Run Python tests and workflow tests](./testing.md#run-python-tests-and-workflow-tests).
4. [Run vignette regression tests](./testing.md#run-vignette-regression-tests).
5. Make fixes to resolve any test failures. If unsure then ask another member of the team for advice.
6. When all tests are passing for your branch, create a new pull request.

Decisions to merge development branches into `develop` are made by Edward Wallace, Felicity Anderson, Kostas Kavoussanakis and Mike Jackon.

## Merging `develop` into `master` and tagging releases

Decisions to merge `develop` into `master` and to tag releases are made by Edward Wallace, Felicity Anderson, Kostas Kavoussanakis and Mike Jackon.
