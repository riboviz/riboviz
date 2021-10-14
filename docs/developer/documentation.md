# Writing and updating documentation

* [Creating Sphinx documentation from Python comments](#creating-sphinx-documentation-from-python-comments)
* [Creating template Sphinx documentation files](#creating-template-sphinx-documentation-files)
* [Updating workflow images](#updating-workflow-images)

---

## Creating Sphinx documentation from Python comments

Create Sphinx pages to reference source code:

```console
$ sphinx-apidoc -f -o py-docs/ riboviz
```

Create HTML documentation

```console
$ cd py-docs
$ make html
```

Open `py-docs/_build/html/index.html` in a browser.

---

## Creating template Sphinx documentation files

The template Sphinx documentation files were originally created as follows:

```console
$ sphinx-quickstart py-docs
> Separate source and build directories (y/n) [n]: y
> Project name: riboviz
> Author name(s): The University of Edinburgh; Rutgers University; University of California, Berkeley
> Project release []: 
> Project language [en]: 
```

Edit `py-docs/conf.py`:

* Uncomment:

```python
# import os
# import sys
```

* Replace:

```python
# sys.path.insert(0, os.path.abspath('.'))
```

* with:

```python
sys.path.insert(0, os.path.abspath('..'))
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

Edit `py-docs/index.rst` and replace content with:

```
riboviz code documentation
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

## Updating workflow images

Workflow images in `docs/images/` are written in the open source [dot](https://graphviz.org/doc/info/lang.html) language from [GraphViz](https://www.graphviz.org/). For an overview, see [Drawing graphs with dot](https://www.graphviz.org/pdf/dotguide.pdf).

If you update the `.dot` files, you should also update the corresponding `.svg` images. GraphViz includes a command-line tool, `dot`, for converting dot files into image in various formats.

To convert a `.dot` file to an `.svg` file:

```console
$ dot -Tsvg workflow.dot > workflow.svg
```

Alternatively, to convert all `.dot` files in `docs/images/`, use the `Makefile`:

```console
$ make clean svgs
```

Alternatively, use an [Editor supporting live preview of GraphViz images](./install.md#editor-supporting-live-preview-of-graphviz-images-optional) that also allows these to be exported as `.svg` images.
