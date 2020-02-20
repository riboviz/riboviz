# Developer guide

## Branching model

The `master` branch is a stable branch. Releases are created via tags from the master branch.

The `develop` branch is for new functionality. At regular intervals this will be merged into the `master` branch.

For branches relating to issues (i.e. new features, significant refactorings or bug fixes), the naming scheme `<name>-<#>` is used, where:

* `<name>` is a short name for the issue. This should be lower-case, with `-` used as a delimiter if desired.
* `<#>` is the number of the issue as in the GitHub issue tracker.

For example, for the issue "Investigate cutadapt -j flag #43" the branch name was `configure-cutadapt-cores-43`.

To request that a branch be merged into the `develop` branch create a new pull request.

---

## Repository structure

```
data/            # Data files used by scripts and vignette
docs/            # Documentation
install/         # Bash scripts to install dependencies
riboviz/         # Python package code
  tools/         # End-user scripts, including prep_riboviz.py
  test/          # pytest-compliant tests
    regression/  # prep_riboviz.py and vignette regression test
rmarkdown/       # Rmarkdown scripts for data preprocessing
rscripts/        # R scripts invoked by vignette
vignette/        # Vignette configuration and input data
website/         # RiboViz Shiny server code and data
```
