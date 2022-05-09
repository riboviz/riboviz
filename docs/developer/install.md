# Install developer dependencies

## Install Python packages

| Package | conda channel | Links |
| ------- | ------------- | ----- |
| pycodestyle | default | [readthedocs](https://pycodestyle.readthedocs.io/), [GitHub](https://github.com/pycqa/pycodestyle) |
| pylint | default | [Pylint](https://www.pylint.org/), [GitHub](https://github.com/PyCQA/pylint/) |
| pytest-cov | default | [pytest-cov](https://pytest-cov.readthedocs.io), [GitHub](https://github.com/pytest-dev/pytest-cov) |
| Sphinx | default | [Sphinx](https://www.sphinx-doc.org/) |

Install:

```console
$ conda install -y pycodestyle
$ conda install -y pylint
$ conda install -y pytest-cov
$ conda install -y sphinx
```

---

## Install R packages

| Package | Links |
| ------- | ----- |
| devtools | [devtools](https://cran.r-project.org/web/packages/devtools/index.html) |
| glue | [glue package on CRAN](https://cran.r-project.org/web/packages/glue/index.html) |
| lintr | [lintr package on CRAN](https://cran.r-project.org/package=lintr), [GitHub](https://github.com/jimhester/lintr) |
| roxygen2 | [roxygen2](https://cloud.r-project.org/web/packages/roxygen2/index.html) |
| styleR | [StyleR package documentation](https://styler.r-lib.org/), [GitHub](https://github.com/r-lib/styler) |
| testthat | [testthat package on CRAN](https://cran.r-project.org/web/packages/testthat/index.html), [testthat package documentation](https://testthat.r-lib.org/) |
| withr | [withr package on CRAN](https://cran.r-project.org/web/packages/withr/index.html) |

Install:

```console
$ R
```
```R
> install.packages("devtools")
> install.packages("glue")
> install.packages("lintr")
> install.packages("roxygen2")
> install.packages("styler")
> install.packages("testthat")
> install.packages("withr")
```

To load the packages before use:

```R
> library(devtools)
> library(glue)
> library(lintr)
> library(roxygen2)
> library(styler)
> library(testthat)
> library(withr)
```

---

## Install `riboviz-py` for development (optional)

To run the Nextflow workflow requires the Python scripts to have been installed via `pip install .`. This installs a `riboviz-py` package, with all riboviz's Python source code into the current Python environment.

If one is working on Python scripts invoked by the workflow it can be time-consuming to repeatedly run `pip install .` after updating the scripts, before rerunning Nextflow.

`pip` allows `riboviz-py` to be installed in "development mode". Rather than installing `riboviz-py`, a link from the Python environment is created back to where you have your Git repository. This means you can edit source files and their changes are immediately available to any code using this Python environment.

`riboviz-py` can be installed in "development mode" as follows:

```console
$ pip install -e .
$ pip list | grep riboviz
riboviz-py                    2.1                 /home/ubuntu/riboviz
```

If you wish to uninstall the "development mode" package in future, run:

```console
$ pip uninstall -y riboviz_py
```

---

## Editor supporting live preview of GraphViz images (optional)

The following free editors supporting live preview of GraphViz images when editing dot documents:

* [Graphviz Support](https://marketplace.visualstudio.com/items?itemName=joaompinto.vscode-graphviz) extension for Microsoft [Visual Studio Code](https://code.visualstudio.com/)
* [GraphViz preview+](https://atom.io/packages/graphviz-preview-plus) for GitHub's [Atom](https://atom.io/).
# Python development
