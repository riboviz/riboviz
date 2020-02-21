# Install developer prerequisites

## Install Python packages

### pylint

Web sites:

* [Pylint](https://www.pylint.org/)
* [BitBucket](https://bitbucket.org/logilab/pylint.org)

Install:

```console
$ conda install -y pylint
```

### pycodestyle

Web sites:

* [readthedocs](https://pycodestyle.readthedocs.io/)
* [GitHub](https://github.com/pycqa/pycodestyle)

Install:

```console
$ conda install -y pycodestyle
```

### pandas

Web sites:

* [pandas](https://pandas.pydata.org/)
* [GitHub](https://github.com/pandas-dev/pandas)

Install:

```console
$ conda install -y pandas
```

### pytest

Web sites:

* [pytest](https://pytest.org/)
* [GitHub](https://github.com/pytest-dev/pytest/)

Install:

```console
$ conda install -y pytest
```

### pytest-cov

Web sites:

* [pytest-cov](https://pytest-cov.readthedocs.io)
* [GitHub](https://github.com/pytest-dev/pytest-cov)

Install:

```console
$ conda install -y pytest-cov
```

### Sphinx

Web site:

* [Sphinx](https://www.sphinx-doc.org/)

```console
$ pip install sphinx
```

---

## Install R packages

### lintr

Web sites:

* [lintr package on CRAN](https://cran.r-project.org/package=lintr)
* [GitHub](https://github.com/jimhester/lintr)

Install:

```console
$ R
```
```R
> install.packages("lintr")
# load the package before use with:
# library("lintr")
```

### styleR

Web sites:

* [StyleR package documentation](https://styler.r-lib.org/)
* [GitHub](https://github.com/r-lib/styler)

Install:

```console
$ R
```
```R
> install.packages("styler")
# load the package before use with:
# library("styler")
```

---

## Install general packages

### GraphViz

Web site: [GraphViz](https://www.graphviz.org/)

Install GraphViz:

**Ubuntu**

```console
$ sudo apt-get install graphviz
```

**CentOS**

```console
$ sudo apt-get install graphviz
```

Check install:

```console
$ dot -V
dot - graphviz version 2.40.1 (20161225.0304)
```

### Editor supporting live preview of GraphViz images (optional)

The following free editors supporting live preview of GraphViz images when editing dot documents:

* [Graphviz Support](https://marketplace.visualstudio.com/items?itemName=joaompinto.vscode-graphviz) extension for Microsoft [Visual Studio Code](https://code.visualstudio.com/)
* [GraphViz preview+](https://atom.io/packages/graphviz-preview-plus) for GitHub's [Atom](https://atom.io/).
