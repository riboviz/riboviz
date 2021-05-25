# Adding and updating dependencies

## Adding a new dependency

If you add a new dependency - a new operating system package, Python package or R package, or another package - then made the following updates.

### Operating system packages

Update [Install general packages](../user/install.md#install-general-packages) in [Install RiboViz and dependencies](../user/install.md):

* Add an entry for the package to the 'Package Links' table.
* Add commands to the 'Install on Ubuntu' and 'Install on CentOS' subsections.
* Add commands to the 'Check packages have installed' subsection.

Add commands to "quick install scripts" (see current content of scripts for what is expected):

* [install-ubuntu.sh](../../bash/install-ubuntu.sh).
* [install-centos.sh](../../bash/install-centos.sh).
* [environment.sh](../../bash/environment.sh), commands to print versions of any command-line tools.
* [environment-tables.sh](../../bash/environment-tables.sh), commands to capture versions of any command-line tools and echo these in a tabular format.

Finally, [Update Dependencies overview](#update-dependencies-overview).

### Bioinformatics tools not available as operating system packages

Update [Install Bioinformatics tools](../user/install.md#install-bioinformatics-tools) in [Install RiboViz and dependencies](../user/install.md):

* Add an entry for the package to the 'Package Links' table.
* Add a new subsection describing how to download, build, configure, install and check the package.
* If the new tool needs the `PATH` configured or other environment variables set, then update commands in `setenv.sh` in [Create `setenv.sh` to configure paths](../user/install.md#create-setenvsh-to-configure-paths).

Add commands to "quick install scripts" (see current content of scripts for what is expected):

* [install-bio-tools.sh](../../bash/install-bio-tools.sh), also update commands to create `setenv.sh` if applicable.
* [environment.sh](../../bash/environment.sh), commands to print versions of any command-line tools.
* [environment-tables.sh](../../bash/environment-tables.sh), commands to capture versions of any command-line tools and echo these in a tabular format.

Finally, [Update Dependencies overview](#update-dependencies-overview).

### Python packages

Update [Install Python packages](../user/install.md#install-python-packages) in [Install RiboViz and dependencies](../user/install.md):

* Add an entry for the package to the 'Package Links' table.
* Add commands to the 'Install' commands.
* Add commands to the 'Check packages have installed command-line tools' commands, if applicable.

Add commands to "quick install scripts" (see current content of scripts for what is expected):

* [install-py.sh](../../bash/install-py.sh).
* [environment-tables.sh](../../bash/environment-tables.sh), add `conda` package names to `CONDA_LIST` or `pip` package names to `PIP_LIST`.

Finally, [Update Dependencies overview](#update-dependencies-overview).

### R packages

Update [Install R and packages required by R packages to be installed](../user/install.md#install-r-and-packages-required-by-r-packages-to-be-installed) in [Install RiboViz and dependencies](../user/install.md):

* Add commands to install any operating system packages to 'Install on Ubuntu' and 'Install on CentOS' subsections.

Update [Install R packages](../user/install.md#install-r-packages) in [Install RiboViz and dependencies](../user/install.md):

* Add an entry for the package to the 'Package Links' table.
* Add commands to the 'Install in R' commands.
* For Bioconductor sub-packages, add commands to the 'For example, for R 3.5 or R 3.6, install in R' and 'For example, for R 3.4, install in R' commands.

Add commands to "quick install scripts" (see current content of scripts for what is expected):

* [install-r-ubuntu.sh](../../bash/install-r-ubuntu.sh), commands to install Ubuntu packages.
 [install-r-centos.sh](../../bash/install-r-centos.sh), commands to install CentOS packages.
* [install-r.R](../../bash/install-r.R), commands to install R packages.
* [install-r-3.4-bioconductor.R](../../bash/install-r-3.4-bioconductor.R), commands to install Bioconductor sub-packages under R 3.4.
* [install-r-3.5-bioconductor.R](../../bash/install-r-3.5-bioconductor.R), commands to install Bioconductor sub-packages under R 3.5+.
* [environment-tables.sh](../../bash/environment-tables.sh), add package names `R_LIST`.

Finally, [Update Dependencies overview](#update-dependencies-overview).

---

## Update an existing dependency

Run through the relevant subsection of [Adding a new dependency](#adding-a-new-dependency), updating the relevant content if necessary.

Finally, [Update Dependencies overview](#update-dependencies-overview).

---

## Update Dependencies overview

To update the [Dependencies overview](#update-dependencies-overview)  tables in [Install RiboViz and dependencies](../user/install.md):

* Run:

```console
$ source bash/environment-tables.sh                                       
```

* Copy the Markdown tables output and paste into `docs/user/install.md`.
