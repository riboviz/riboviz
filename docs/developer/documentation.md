# Writing and updating documentation

## Create Sphinx documentation from Python

Create Sphinx pages to reference source code:

```console
$ sphinx-apidoc -f -o code-docs/ riboviz
```

Create HTML documentation

```console
$ cd code-docs
$ make html
```

Open `code-docs/_build/html/index.html` in a browser.

---

## How the template Sphinx documentation files were originally created

```console
$ sphinx-quickstart code-docs
> Separate source and build directories (y/n) [n]: y
> Project name: RiboViz
> Author name(s): The University of Edinburgh; Rutgers University; University of California, Berkeley
> Project release []: 
> Project language [en]: 
```

Edit `code-docs/conf.py`:

* Uncomment:

```python
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
```

* Update above to:

```python
sys.path.insert(0, os.path.abspath('../../riboviz'))
```

* Update:

```python
extensions = [
]
```

* to:

```python
extensions = [
    'sphinx.ext.autodoc'
]
```

Edit `code-docs/index.rst` and replace content with:

```
RiboViz code documentation
==========================

.. toctree::
   :maxdepth: 1
   :caption: Code documentation

   modules

Indices and tables:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```

---

## Update workflow images

Workflow images in [images](../images) are written in the open source [dot](https://graphviz.org/doc/info/lang.html) language from [GraphViz](https://www.graphviz.org/). For an overview, see [Drawing graphs with dot](https://www.graphviz.org/pdf/dotguide.pdf). GraphViz includes a command-line tool, `dot`, for converting dot files into image in various formats.

Convert `.dot` file to `.svg` file:

```console
$ dot -Tsvg workflow.dot > workflow.svg
```

When `.dot` files are updated, the corresponding `.svg` images should be updated also.
