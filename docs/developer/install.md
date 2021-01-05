# Install developer dependencies

## Install Python packages

| Package | Links |
| ------- | ----- |
| pylint | [Pylint](https://www.pylint.org/), [BitBucket](https://bitbucket.org/logilab/pylint.org) |
| pycodestyle | [readthedocs](https://pycodestyle.readthedocs.io/), [GitHub](https://github.com/pycqa/pycodestyle) |
| pytest-cov | [pytest-cov](https://pytest-cov.readthedocs.io), [GitHub](https://github.com/pytest-dev/pytest-cov) |
| Sphinx | [Sphinx](https://www.sphinx-doc.org/) |

Install:

```console
$ conda install -y pylint
$ conda install -y pycodestyle
$ conda install -y pytest-cov
$ pip install sphinx
```

---

## Install R packages

| Package | Links |
| ------- | ----- |
| lintr | [lintr package on CRAN](https://cran.r-project.org/package=lintr), [GitHub](https://github.com/jimhester/lintr) |
| styleR | [StyleR package documentation](https://styler.r-lib.org/), [GitHub](https://github.com/r-lib/styler) |
| roxygen2 | [roxygen2](https://cloud.r-project.org/web/packages/roxygen2/index.html) |


```console
$ R
```
```R
> install.packages("lintr")
> install.packages("styler")
> install.packages("roxygen2")
```

To load the packages before use:

```R
> library(lintr)
> library(styler)
> library(roxygen2)
```

---

## Install general packages

| Package | Links |
| ------- | ----- |
| GraphViz | [GraphViz](https://www.graphviz.org/) |

Install:

**Ubuntu**

```console
$ sudo apt-get install graphviz
```

**CentOS**

```console
$ sudo yum install graphviz
```

Check install:

```console
$ dot -V
```

---

## Editor supporting live preview of GraphViz images (optional)

The following free editors supporting live preview of GraphViz images when editing dot documents:

* [Graphviz Support](https://marketplace.visualstudio.com/items?itemName=joaompinto.vscode-graphviz) extension for Microsoft [Visual Studio Code](https://code.visualstudio.com/)
* [GraphViz preview+](https://atom.io/packages/graphviz-preview-plus) for GitHub's [Atom](https://atom.io/).
