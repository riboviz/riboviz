# Install developer dependencies

## Install Python packages

| Package | conda channel | Links |
| ------- | ------------- | ----- |
| pycodestyle | default | [readthedocs](https://pycodestyle.readthedocs.io/), [GitHub](https://github.com/pycqa/pycodestyle) |
| pylint | default | [Pylint](https://www.pylint.org/), [BitBucket](https://bitbucket.org/logilab/pylint.org) |
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

## Editor supporting live preview of GraphViz images (optional)

The following free editors supporting live preview of GraphViz images when editing dot documents:

* [Graphviz Support](https://marketplace.visualstudio.com/items?itemName=joaompinto.vscode-graphviz) extension for Microsoft [Visual Studio Code](https://code.visualstudio.com/)
* [GraphViz preview+](https://atom.io/packages/graphviz-preview-plus) for GitHub's [Atom](https://atom.io/).
