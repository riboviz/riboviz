# Adding and updating dependencies

* [Adding new operating system packages](#adding-new-operating-system-packages)
* [Adding tools not available as operating system packages](#adding-tools-not-available-as-operating-system-packages)
* [Adding Python packages](#adding-python-packages)
* [Adding R packages](#adding-r-packages)
* [Updating existing dependencies](#updating-existing-dependencies)
* [Updating Dependencies overview tables](#updating-dependencies-overview-tables)

---

## Adding new operating system packages

Update 'Install operating system packages' in `docs/user/install.md`:

* Add an entry for the package to the 'Package Links' table.
* Add commands to the 'Install on Ubuntu' and 'Install on CentOS' subsections.
* Add commands to the 'Check packages have installed' subsection.

Add commands to "quick install scripts" (see current content of scripts for what is expected):

* `bash/install-ubuntu.sh`
* `bash/install-centos.sh`
* `bash/environment-tables.sh`, add commands to capture versions of any command-line tools and echo these in a tabular format.

Finally, see [Updating Dependencies overview tables](#updating-dependencies-overview-tables).

---

## Adding tools not available as operating system packages

Update 'Install tools not available as operating system packages' in `docs/user/install.md`:

* Add an entry for the package to the 'Package Links' table.
* Add a new subsection describing how to download, build, configure, install and check the package.
* If the new tool needs the `PATH` configured or other environment variables set, then update commands in `setenv.sh` in 'Create `setenv.sh` to configure paths'.

Add commands to "quick install scripts" (see current content of scripts for what is expected):

* `bash/install-tools.sh`, also update commands to create `setenv.sh`, if applicable.
* `bash/environment-tables.sh`, add commands to capture versions of any command-line tools and echo these in a tabular format.

Finally, [Updating Dependencies overview tables](#updating-dependencies-overview-tables).

---

## Adding Python packages

Update 'Install Python packages' in `docs/user/install.md` or `docs/developer/install.md` (depending on whether the package is for all users or for developers only):

* Add an entry for the package to the 'Package Links' table.
* Add commands to the 'Install' commands.
* For user-specific packages, add commands to the 'Check packages have installed command-line tools' commands, if applicable.

Add commands to "quick install scripts" (see current content of scripts for what is expected):

* `bash/install-py.sh`
* `bash/environment-tables.sh`, add `conda` package names to `CONDA_LIST` or `pip` package names to `PIP_LIST`.

Finally, [Updating Dependencies overview tables](#updating-dependencies-overview-tables).

---

## Adding R packages

Update 'Install R and packages required by R packages to be installed' in `docs/user/install.md`:

* Add commands to install any operating system packages to 'Install on Ubuntu' and 'Install on CentOS' subsections.

Update 'Install R packages' in `docs/user/install.md` or `docs/developer/install.md` (depending on whether the package is for all users or for developers only):

* Add an entry for the package to the 'Package Links' table.
* Add commands to the 'Install in R' commands.
* For Bioconductor sub-packages, add commands to the 'For example, for R 3.5 or R 3.6, install in R' and 'For example, for R 3.4, install in R' commands.

Add commands to "quick install scripts" (see current content of scripts for what is expected):

* `bash/install-r-ubuntu.sh`, commands to install Ubuntu packages.
* `bash/install-r-centos.sh`, commands to install CentOS packages.
* `bash/install-r.R`, commands to install R packages.
* `bash/install-r-3.4-bioconductor.R`, commands to install Bioconductor sub-packages under R 3.4.
* `bash/install-r-3.5-bioconductor.R`, commands to install Bioconductor sub-packages under R 3.5+.
* `bash/environment-tables.sh`, add package names `R_LIST`.

Finally, [Updating Dependencies overview tables](#updating-dependencies-overview-tables).

---

## Updating existing dependencies

Run through the section corresponding to the nature of the dependency being updated and make any required updates to documentation and scripts:

* [Adding new operating system packages](#adding-new-operating-system-packages)
* [Adding tools not available as operating system packages](#adding-tools-not-available-as-operating-system-packages)
* [Adding Python packages](#adding-python-packages)
* [Adding R packages](#adding-r-packages)

Finally, [Updating Dependencies overview tables](#updating-dependencies-overview-tables).

---

## Updating Dependencies overview tables

To update the 'Dependencies overview' tables in `docs/user/install.md`:

* Run:

```console
$ source bash/environment-tables.sh                                       
```

* Copy the Markdown tables that are output and paste into `docs/user/install.md`.
