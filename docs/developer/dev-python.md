# Developing Python components

* [Python style](#python-style)
  - [File names](#file-names)
  - [Comments](#comments)
  - [Style checking Python code](#style-checking-python-code)
* [Python command-line tools](#python-command-line-tools)
  - [Defining command-line parameters](#defining-command-line-parameters)
  - [Adding, renaming, and removing command-line tools](#adding-renaming-and-removing-command-line-tools)
* [Comments and documentation](#comments-and-documentation)
  - [Creating Sphinx documentation from Python comments](#creating-sphinx-documentation-from-python-comments)
  - [Creating template Sphinx documentation files](#creating-template-sphinx-documentation-files)
* [Testing](#testing)
  - [Run Python tests and workflow tests](#run-python-tests-and-workflow-tests)
  - [Useful pytest flags](#useful-pytest-flags)
  - [Information on pytest fixtures and parameters](#information-on-pytest-fixtures-and-parameters)
* [Miscellaneous](#miscellaneous)
  - [Clone conda environments](#clone-conda-environments)
  - [Understanding YAML `NULL` and Python `None`](#understanding-yaml-null-and-python-none)

---

## Python style

Python code should be formatted to conform as far as possible to the [PEP 8 -- Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/).

The tools [pylint](https://pylint.org/) and [pycodestyle](https://pycodestyle.pycqa.org/en/latest/) can help by highlighting deviations from the coding style. See [Style checking Python code](#style-checking-python-code) below.

See also the riboviz [Style guide](./style-guide.md).

### File names

Python file names must be in snake-case i.e., lower-case and delimited by underscores, not hyphens. For example, `sample_sheets.py`, `get_cds_codons.py`.

### Comments

Python code should be commented using [doc-strings](https://www.python.org/dev/peps/pep-0257/). The doc-strings should include [Sphinx](https://www.sphinx-doc.org/en/master/)-compliant [reStructuredText](https://docutils.sourceforge.io/rst.html) markup. Examples are given below.

Modules should have a doc-string descring the module. For example:

```python
"""
Workflow configuration parameter names.
"""
```

Constants and module-level variable should have doc-strings briefly describing what they are. For example:

```python
INPUT_DIR = "dir_in"
""" Input directory. """
```

Functions should have doc-strings describing what they do and their parameters, return values (if applicable) and any exceptions (if applicable). Types of parameters and return values should also be documented. For example:

```python
def load_sample_sheet(file_name, delimiter="\t", comment="#"):
    """
    Load a sample sheet from a file. The sample sheet is assumed to
    have a header with column names ``SampleID`` and ``TagRead``.

    :param file_name: File name
    :type file_name: str or unicode
    :param delimiter: Delimiter
    :type delimiter: str or unicode
    :param comment: Comment prefix
    :type comment: str or unicode
    :return: Sample sheet
    :rtype: pandas.core.frame.DataFrame
    :raise FileNotFoundError: If the file cannot be found or is \
    not a file
    :raise AssertionError: If there is no header with ``SampleID`` \
    and ``TagRead`` columns.
    """
```

To use fixed-width (teletype-style) text use ` `` ` (see the markup on `SampleID` and `TagRead` above).

For references to other Python modules, functions or constants use Sphinx's `:py:` roles. For example:

* Cross-reference to a module (`:py:mod:`):

```python
Count reads using :py:mod:`riboviz.tools.count_reads`.
```

* Cross-reference to a constant (`:py:constant:`):

```python
... also matching :py:const:`riboviz.workflow_files.UMI_EXTRACT_FQ` ...
```

* Cross-reference to a function (`:py:func:`):

```python
See :py:func:`riboviz.count_reads.count_reads` for information ...
```

Note that if the cross-reference is to an entity in the same module then only the local name of the entity needs to be specified. For example:

```python
See :py:func:`count_reads` for information ...
```

For examples and further information, see:

* Sphinx's Python [Info field lists](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#info-field-lists) on the Python fields that Sphinx can recognise.
* Sphinx's [Cross-referencing Python objects](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#cross-referencing-python-objects) on the available cross-references.
* Sphinx's [The Python Domain](https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#the-python-domain).

See below for information on [Creating Sphinx documentation from Python comments](#creating-sphinx-documentation-from-python-comments).

### Style checking Python code

[pylint](https://pylint.org/) and [pycodestyle](https://pycodestyle.pycqa.org/en/latest/) can be run on individual files, groups of files or every file in a directory and its subdirectories. For example:

```console
$ pylint riboviz/check_fasta_gff.py
$ pycodestyle riboviz/check_fasta_gff.py
```
```console
$ pylint riboviz/
$ pycodestyle riboviz/
```

---

## Python command-line tools

Python command-line tools must be placed within `riboviz/tools/` i.e. within the `riboviz.tools` module.

Command-line parsing must be implemented using [argparse](https://docs.python.org/3/library/argparse.html).

`riboviz.tools` scripts must not contain any functionality beyond defining and parsing command-line parameters then invoking the functionality of the script. The functionality of the script itself should be placed in module(s) under `riboviz`.

`riboviz.tools` scripts must contain an entry point defined as follows:

```
if __name__ == "__main__":
    main()
```

where the `main()` function is a function in the script responsible for parsingthe command-line parameters and invoking the functionality of the script itself, as defined in other `riboviz` module(s). See the existing `riboviz.tools` implementations for the pattern to use.

### Defining command-line parameters

The `dest` parameter of `ArgumentParser.add_argument` can be use to explicitly define the Python variables into which a command-line parameter is to be placed. For example:

```python
parser = argparse.ArgumentParser(description="Program")
parser.add_argument("-o", "--output-dir", dest="output_dir", nargs='?',
                    help="Output directory")
options = parser.parse_args()
output_dir = options.output_dir
```

### Adding, renaming, and removing command-line tools

If adding, renaming or removing	command-line tools:

* Update 'Python command-line tools' table in `docs/user/command-line-tools.md`.
* Update `console_scripts` value  in `[options.entry_points]` in `setup.cfg`.

---

## Comments and documentation

### Creating Sphinx documentation from Python comments

Create Sphinx pages to reference source code:

```console
$ sphinx-apidoc -f -o py-docs/ riboviz
```

Create HTML documentation

```console
$ cd py-docs
$ make html
```

If there are errors in the comment formatting then these will be displayed.

Open `py-docs/_build/html/index.html` in a browser.

Manually check the rendered comments to see that they have no errors that cannot be caught during the build. These include whether or not you have:

* Marked up ` ``teletype font`` ` text correctly.
* Defined cross-references to modules, constants, and functions correctly.

### Creating template Sphinx documentation files

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

## Testing

### Run Python tests and workflow tests

Run:

```console
$ pytest --ignore-glob="*integration*"
```

**Troubleshooting: `PendingDeprecationWarning`**

`PendingDeprecationWarning` `warnings` can be ignored.

### Useful pytest flags

* `-s`: disable output capture so, for example, `print` messages are shown.
* `-v`: verbose mode, displays names of test functions run.
* `-k`: run a specific test function.
* `--cov-config=.coveragerc --cov-report term-missing --cov=riboviz`: create a test coverage report which includes the line numbers of statements that were not executed.

### Information on pytest fixtures and parameters

The riboviz integration and Python unit tests make extensive use of pytest fixtures and parameterised tests. For more information on pytest fixtures and parameters see:

* [pytest fixtures: explicit, modular, scalable](https://docs.pytest.org/en/6.2.x/fixture.html)
* [Parametrizing fixtures and test functions](https://docs.pytest.org/en/6.2.x/parametrize.html)
* [Basic patterns and examples](https://docs.pytest.org/en/6.2.x/example/simple.html)

---

## Miscellaneous

### Clone conda environments

If you want to install a new Python package to explore it without changing your current conda environment you can first clone your conda environment, then install the package into the clone. For example:

```console
$ conda create --name riboviz-test-install --clone riboviz
$ conda activate riboviz-test-install
```

To remove an environment when you are finished with it you can run, for example:

```console
$ conda deactivate riboviz-test-install
$ conda env remove --name riboviz-test-install
```

For more information, see conda's [Managing environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

### Understanding YAML `NULL` and Python `None`

If a parameter in a YAML file has value `null`, `NULL` or no value at all then, after reading the file into Python (using the `yaml` library), it will have value `None`. For example, given a YAML file with:

```yaml
a: null
b: NULL
c:
```

Loading this with

```python
with open("config.yml") as f:
    config = yaml.load(f, Loader=yaml.SafeLoader)
```

would result in the following all having value `None`:

```python
config['a']
config['b']
config['c']
```
