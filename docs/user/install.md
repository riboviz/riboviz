# Install RiboViz and dependencies

## About these instructions

These instructions written for Ubuntu and CentOS and tested upon Ubuntu 18.04 and CentOS 7. Other Linux flavours will require different commands to be run.

You need to have permission to run `sudo` to install many of the dependencies. If you don't have `sudo` access you will have to ask a local system administrator to run these commands for you.

### Windows users

We suggest that you:

* Either, use a virtual machine running under [VMWare Workstation Player](https://www.vmware.com/uk/products/workstation-player.html) or [Oracle VirtualBox](https://www.virtualbox.org/). We provide quick-start instructions for using VMWare to:
  - [Deploy a Ubuntu Virtual Machine using VMWare on Windows 10](./deploy-ubuntu-vmware-windows.md).
  - [Deploy a CentOS Virtual Machine using VMWare on Windows 10](./deploy-centos-vmware-windows.md).
* Or, use Windows Subsystem for Linux following [Windows Subsystem for Linux Installation Guide for Windows 10](https://docs.microsoft.com/en-us/windows/wsl/install-win10). We have successfully installed and run RiboViz under WSL2 using Microsoft Store's [Ubuntu 18.04 LTS](https://www.microsoft.com/en-gb/p/ubuntu-1804-lts/9n9tngvndl3q).

### Mac OS users

We suggest that you check out the web sites for each prerequisite for information on how to install the prerequisites under Mac OS.

---

## Dependencies overview

The following tables summarise the packages required by RiboViz. Instructions to install each dependency are given in the following sections. Only minimal installation instructions are given. For full information, see the documentation for each dependency.

The versions listed are those used by a RiboViz developer when preparing the current release. Other versions may also be acceptable, but see the constraints below.

| Command-line tool | Version |
| ----------------- | ------- |
| Git | 2.17.1 |
| cURL | 7.69.1 |
| bedtools | 2.26.0 |
| hdf5tools (h5diff) | 1.10.6 |
| pigz | 2.4 |
| pandoc | 1.19.2.4 |
| GraphViz (dot) | 2.40.1 |
| zip | 3.0 |
| unzip | 6.00 |
| R | 3.6.3 |
| Python | 3.7.6 |
| Cutadapt | 1.18 |
| samtools | 1.9 |
| UMI-tools | 1.0.1 |
| Java (javac) | 1.8.0_152-release |
| Java (java) | 1.8.0_152-release |
| Nextflow | 20.04.1.5335 |
| Hisat2 | 2.1.0 |
| Bowtie | 1.2.2 |

| R Package | Version |
| --------- | ------- |
| Biostrings | 2.54.0 |
| devtools | 2.3.2 |
| git2r | 0.27.1 |
| glue | 1.4.1 |
| here | 1.0.1 |
| knitr | 1.33 |
| lintr | 2.0.1 |
| optparse | 1.6.6 |
| plotly | 4.9.2.1 |
| RcppRoll | 0.3.0 |
| rhdf5 | 2.30.1 |
| rmarkdown | 2.7 |
| roxygen2 | 7.1.1 |
| Rsamtools | 2.2.3 |
| rtracklayer | 1.46.0 |
| shiny | 1.5.0 |
| ShortRead | 1.44.3 |
| styler | 1.3.2 |
| testthat | 3.0.1 |
| tidyverse | 1.3.0 |
| withr | 2.3.0 |

| Python Package | Version | Package Manager |
| -------------- | ------- | --------------- |
| biopython | 1.77 | conda |
| cutadapt | 1.18 | conda |
| gffutils | 0.10.1 | conda |
| gitpython | 3.1.3 | conda |
| h5py | 2.10.0 | conda |
| nextflow | 20.04.1 | conda |
| pandas | 1.0.5 | conda |
| pycodestyle | 2.6.0 | conda |
| pylint | 2.5.3 | conda |
| pysam | 0.15.3 | conda |
| pytest | 5.4.3 | conda |
| pytest-cov | 2.10.1 | conda |
| pyyaml | 5.3.1 | conda |
| samtools | 1.9 | conda |
| umi_tools | 1.0.1 | conda |
| sphinx | 3.2.0 | pip |

Certain packages are only required if you plan to develop and extend RiboViz. These packages are (see [Install developer dependencies](../developer/install.md)):

* R: devtools, glue, lintr, roxygen2, styler, testthat, withr.
* Python: pycodestyle, pylint, pytest-cov, sphinx.

Requirements and constraints:

* R 3.6, or later, is required.
* Python 3 is required. Python 2 came to the end of its supported life in 2020 and there will be no Python 2.8 (see [PEP 373 Python 2.7 Release Schedule](https://legacy.python.org/dev/peps/pep-0373/)).
* Python 3.6, or later, is strongly recommended.
* Either [Miniconda](https://conda.io/miniconda.html) or [Anaconda Distribution](https://www.anaconda.com/distribution/) are strongly recommended.
* Cutadapt 1.18 (2018-09-07), or later, is required.
* Samtools 1.9 is required if running on CentOS 7.
* Nextflow 20, or later, is required.
* Hisat 2.1.0 is strongly recommended. Hisat2 2.2.0 users have reported bugs and issues (see for example [DaehwanKimLab/hisat2#242](https://github.com/DaehwanKimLab/hisat2/issues/242) and [DaehwanKimLab/hisat2#245](https://github.com/DaehwanKimLab/hisat2/issues/245)) which Hisat2 will resolve in a future release.

---

## Install operating system packages

| Package | Links |
| ------- | ----- |
| Git | [Git](https://git-scm.com/) |
| cURL | [cURL](https://curl.haxx.se/) |
| EPEL (CentOS only) | [EPEL](https://fedoraproject.org/wiki/EPEL) |
| bedtools | [bedtools](http://bedtools.readthedocs.io/en/latest/), [GitHub](https://github.com/arq5x/bedtools2) |
| hdf5tools | [HDF5](https://portal.hdfgroup.org/display/HDF5) |
| pigz | [pigz](http://zlib.net/pigz/) |
| pandoc | [pandoc](https://pandoc.org) |
| GraphViz | [GraphViz](https://www.graphviz.org/) |

plus common utilities and low level libraries.

*Note:** It is OK if some of these packages are already present.

### Install on Ubuntu

```console
$ sudo apt update -y
$ sudo apt install -y git
$ sudo apt install -y curl
$ sudo apt install -y bedtools
$ sudo apt install -y hdf5-tools
$ sudo apt install -y pigz
$ sudo apt install -y pandoc
$ sudo apt install -y graphviz
$ sudo apt install -y zip
$ sudo apt install -y unzip
$ sudo apt install -y libxml2-dev
$ sudo apt install -y libssl-dev
$ sudo apt install -y libcurl4-openssl-dev
$ sudo apt install -y libgit2-dev
```

### Install on CentOS

```console
$ sudo yum update -y
$ sudo yum install -y git
$ sudo yum install -y curl
$ sudo yum install -y epel-release
$ sudo yum install -y BEDTools
$ sudo yum install -y hdf5-devel
$ sudo yum install -y pigz
$ sudo yum install -y pandoc
$ sudo yum install -y graphviz
$ sudo yum install -y libxml2-devel
$ sudo yum install -y openssl-devel
$ sudo yum install -y libcurl-devel
$ sudo yum install -y libjpeg-devel
$ sudo yum install -y libgit2-devel
```

### Check tools have installed

```console
$ git --version
$ curl --version
$ bedtools -version
$ h5diff --version
$ pigz --version
$ pandoc --version
$ dot -V
$ zip -v
$ unzip -v
```

`h5diff` is one of the hdf5tools.

---

## Install R 3.6+

| Package | Links |
| ------- | ----- |
| R | [The R Project for Statistical Computing](https://www.r-project.org/), [The Comprehensive R Archive Network](https://cran.r-project.org/) (CRAN) |

**Note:** R 3.6, or later, is required.

### Install R and packages required by R packages to be installed

**Install on Ubuntu**

```console
$ sudo apt install -y r-base
$ sudo apt install -y r-base-dev
```
```console
$ R --version
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
```

Your version of R can differ from that shown but must be version 3.6 or above.

Default package managers may not have the most up-to-date version of R available. This may cause problems when installing R or its packages. [The Comprehensive R Archive Network](https://cran.r-project.org/) has information on alternative ways to get a more recent version of R.

For Ubuntu, you can do the following, from [Ubuntu Packages For R - Brief Instructions](https://cran.r-project.org/bin/linux/ubuntu/):

```console
$ sudo apt install -y --no-install-recommends software-properties-common dirmngr
$ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
$ sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
$ sudo apt update -y
$ sudo apt install -y r-base
$ sudo apt install -y r-base-dev
$ R --version
R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
```

**Install on CentOS**

```console
$ sudo yum install -y R
$ sudo yum install -y R-devel
```
```console
$ R --version
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
```

Your version of R can differ from that shown but must be version 3.6 or above.

If you wish to install the latest version of R then [RStudio](https://docs.rstudio.com/) has instructions on how to do this at [Install R](https://docs.rstudio.com/resources/install-r/). For example, to install R 4.1.0:

```console
$ sudo yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm 
$ sudo yum-config-manager --enable "rhel-*-optional-rpms"
$ export R_VERSION=4.1.0
$ curl -O https://cdn.rstudio.com/r/centos-7/pkgs/R-${R_VERSION}-1-1.x86_64.rpm
$ sudo yum install -y R-${R_VERSION}-1-1.x86_64.rpm
$ /opt/R/${R_VERSION}/bin/R --version
R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
$ sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
$ sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript
$ R --version
R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
```

---

## Install R packages

| Package | Links |
| ------- | ----- |
| RcppRoll | [RcppRoll](https://cran.r-project.org/web/packages/RcppRoll/index.html) |
| git2r | [git2r](https://docs.ropensci.org/git2r), [GitHub](https://github.com/ropensci/git2r) |
| here | [here](https://here.r-lib.org/), [CRAN](https://cran.r-project.org/package=here), [GitHub](https://github.com/r-lib/here) |
| knitr | [knitr](https://cran.r-project.org/web/packages/knitr/index.html) |
| optparse | [optparse](https://cran.r-project.org/web/packages/optparse/index.html) |
| plotly | [plotly](https://cran.r-project.org/web/packages/plotly/index.html) |
| rmarkdown | [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html) |
| shiny | [shiny](https://cran.r-project.org/web/packages/shiny/index.html) |
| tidyverse | [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) |
| Bioconductor Biostrings | [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) |
| Bioconductor Rsamtools | [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) |
| Bioconductor ShortRead | [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) |
| Bioconductor rhdf5 | [rhdf5](https://bioconductor.org/packages/release/bioc/html/rhdf5.html) |
| Bioconductor rtracklayer | [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) |

**Note:** Biostrings is installed as a dependency of Rsamtools.

If, when installing R packages, you see a message like:

```
Would you like to use a personal library instead? (yes/No/cancel) yes
Would you like to create a personal library
"~/R/x86_64-pc-linux-gnu-library/3.6"
to install packages into? (yes/No/cancel)
```

then enter `yes`.

Install in R:

```console
> R
```
```r
> install.packages("RcppRoll")
> install.packages("git2r")
> install.packages("here")
> install.packages("knitr")
> install.packages("optparse")
> install.packages("plotly")
> install.packages("rmarkdown")
> install.packages("shiny")
> install.packages("tidyverse")
```

The commands to install Bioconductor packages depend on your version of R. For full details:

* See Bioconductor/R compatibility on [Bioconductor releases](https://bioconductor.org/about/release-announcements/).
* Click the link of a Bioconductor release consistent with your version of R.
* Click the link of the specific package.

For example, for R 3.6 or R4.1, install in R:

```r
> install.packages("BiocManager")
> BiocManager::install("Rsamtools")
> BiocManager::install("rtracklayer")
> BiocManager::install("rhdf5")
> BiocManager::install("ShortRead")
```

### Troubleshooting: `installation path not writeable`

The following warning can be ignored:

```
installation path not writeable, unable to update packages: foreign
```

See [Question: unable to update packages: foreign, Matrix](https://support.bioconductor.org/p/96834/) for an explanation:

> That's not an error! You just got an informative message saying that two of the base R packages couldn't be updated...
> All of your Bioconductor packages end up in the first dir, which is writeable by you, and the base and core packages go in the second dir, which is only writeable by an administrator. You shouldn't be running R as an administrator, like ever, so it's common for you to get the message that you saw. If you really care to update the core packages, you can run R as an administrator, do biocLite, and then restart as a lower-permissioned user after the update.

You can check that it is available (your exact path may differ):

```console
$ ls ~/R/x86_64-pc-linux-gnu-library/3.4/Rsamtools/
DESCRIPTION  libs  LICENSE  Meta  NAMESPACE  NEWS
```

```console
$ ls ~/R/x86_64-redhat-linux-gnu-library/3.6/Rsamtools/
DESCRIPTION  libs  LICENSE  Meta  NAMESPACE  NEWS
```

### Troubleshooting: `Cannot allocate memory`

```
Error in system2(file.path(R.home("bin"), "R"), c(if (nzchar(arch))
paste0("--arch=",  :
  cannot popen ' '/usr/lib/R/bin/R' --no-save --slave 2>&1 <
  '/tmp/Rtmpw3pOH7/file12471113d0d2b'', probable reason 'Cannot
  allocate memory'
* removing "/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.6/Rsamtools"
Warning in q("no", status = status, runLast = FALSE) :
  system call failed: Cannot allocate memory

The downloaded source packages are in
        "/tmp/RtmpOEVbaL/downloaded_packages"
installation path not writeable, unable to update packages: foreign
Warning message:e
In install.packagees(pkgs = doing, lib = lib, ...) :
  installation of package "Rsamtools" had non-zero exit status
```

You may need to assign more memory to R or your machine.

### Troubleshooting: `package "XML" is not available`

If you get this error message when running:

```R
> BiocManager::install("rtracklayer")
```

then one solution may be to install "XML" specifying the URL of the source package, for example:

```R
> install.packages("https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz", repos=NULL, type="source")
```

### Troubleshooting: `ShortRead` installation fails on CentOS

If the following failure arises when installing `ShortRead` under CentOS:

```R
> install.packages("BiocManager")
...
> BiocManager::install("ShortRead")
Warning messages:
1: In .inet_warning(msg) :
  installation of package "png" had non-zero exit status
2: In .inet_warning(msg) :
  installation of package "jpeg" had non-zero exit status
3: In .inet_warning(msg) :
  installation of package "latticeExtra" had non-zero exit status
4: In .inet_warning(msg) :
  installation of package "ShortRead" had non-zero exit status
Install libjpegdevel:
```

Failure 2 causes failure 4 may be due to a missing package. To resolve this failure, install the `libjpeg-devel` package:

```console
$ sudo yum install -y libjpeg-devel
```

Failure 1 causes failures 3 and 4 and may be due to, on CentOS 7, the system-wide `libpng` package being at version 15, while Miniconda 3, after installation of the Python packages above, has `libpng` version 16. To resolve this failure, start a new bash terminal but do not activate the Miniconda environment, start R and reinstall the package:

```console
$ R
```
```R
> install.packages("BiocManager")
> BiocManager::install("ShortRead")
```

---

## Install Python

| Package | Links |
| ------- | ----- |
| Python | [python](https://www.python.org/) |

The instructions which follow have been written for Miniconda Python. If using Anaconda then, when installing some packages, you will be told that they are already available. This is because Anaconda comes with a wide range of common Python packages.

If you are using other distributions of Python, you will need to consult the relevant documentation for each package for installation information. See also the section on [python and python3](#python-and-python3) below.

**Note:** Python 3 is required. Python 2 came to the end of its supported life in 2020 and there will be no Python 2.8 (see [PEP 373 Python 2.7 Release Schedule](https://legacy.python.org/dev/peps/pep-0373/)).

**Note:** Either [Miniconda](https://conda.io/miniconda.html) or [Anaconda Distribution](https://www.anaconda.com/distribution/) are strongly recommended.

**Note:** Python 3.6, or later, is strongly recommended.

### Install Miniconda Python 3.6+

On Linux, install:

```console
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
$ bash miniconda3.sh -b -p $HOME/miniconda3
```

**Note:** Make sure you use `-O`, which provides a name for the downloaded file, and not `-o`, which provides the name of a file for messages about the download.

For Mac OS and Windows installers, go to [Miniconda installation page](https://docs.conda.io/en/latest/miniconda.html).

When Miniconda has installed, activate the environment:

```console
$ source $HOME/miniconda3/bin/activate
```

Create a `riboviz` environment and activate it:

```console
$ conda create -y --name riboviz python=3.7
$ conda activate riboviz
$ python --version
```

Your version of Python may differ from that shown.

**Troubleshooting: `...command not found...`**

If you see:

```
$ bash miniconda3.sh -b -p $HOME/miniconda3
miniconda3.sh: line 1: --2019-07-31: command not found
miniconda3.sh: line 2: syntax error near unexpected token `('
miniconda3.sh: line 2: `Resolving repo.continuum.io
(repo.continuum.io)... 104.18.200.79, 104.18.201.79,
2606:4700::6812:c94f, ...'
```

then rerun `wget` and use `-O`, not `-o`.

**Troubleshooting: incompatible Python versions**

If you find when installing packages below that your version of Python is too new (e.g. your Python version is 3.8 and the tool only works with Python 3.7), then two options are:

1. See if a more recent version of the tool is available. For example a version available via `pip` may be more up-to-date than a version available via `conda`.
2. Create a conda environment that supports a version of Python compatible with the tools.

```console
$ conda create -y --name riboviz python=<VERSION>
```

  - For example:

```console
$ conda create --name riboviz python=3.7
```

### `python` and `python3`

If you have an environment which has both `python` and `python3`, such as can arise when you are using a system that has both Python 2 and Python 3 packages centrally installed, then please note the following.

The RiboViz workflow invokes both Python and R scripts. It invokes Python scripts using the command `python`. If you have a system that has both `python`, which invokes Python 2, and `python3`, which invokes Python 3, then the workflow will fail as RiboViz's Python scripts are Python 3-compatible only.

A workaround is to create a local `bin` directory with a symbolic link called `python` which points to Python 3 (and similarly for `pip` and `pip3`). This can be done as follows:

```console
$ mkdir ~/bin
$ cd ~/bin
$ ln -s $(which python3) python
$ ln -s $(which pip3) pip
$ cd
```

Now, when you run `python`, `python3` should be invoked. If the symlinks aren't picked up then you may need to add ~/bin to your PATH:

```console
$ export PATH=~/bin:$PATH
```

This approach was suggested in a [comment](https://stackoverflow.com/a/55295939) on StackOverflow's [Unable to set default python version to python3 in ubuntu](https://stackoverflow.com/questions/41986507/unable-to-set-default-python-version-to-python3-in-ubuntu).

We would recommend using Miniconda, Anaconda or some other virtual environment solution for Python which provide a more usable means of managing multiple environments (including Python 2 and Python 3).

---

## Install Python and conda packages

| Package | Links |
| ------- | ----- |
| pyyaml | [PyYAML](https://pyyaml.org/), [GitHub](https://github.com/yaml/pyyaml/) |
| gitpython | [gitpython](https://gitpython.readthedocs.io/en/stable/), [GitHub](https://github.com/gitpython-developers/GitPython) |
| pytest | [pytest](https://pytest.org/), [GitHub](https://github.com/pytest-dev/pytest/) |
| pandas | [pandas](https://pandas.pydata.org/), [GitHub](https://github.com/pandas-dev/pandas) |
| Cutadapt | [GitHub](https://github.com/marcelm/cutadapt), [readthedocs](https://cutadapt.readthedocs.io/) |
| pysam | [GitHub](https://github.com/pysam-developers/pysam/), [readthedocs](https://pysam.readthedocs.io/) |
| Samtools | [samtools](https://www.htslib.org/) |
| BioPython | [Biopython](http://biopython.org/) |
| gffutils | [gffutils](http://daler.github.io/gffutils/) |
| h5py | [h5py](https://www.h5py.org/) |
| UMI-tools | [GitHub](https://github.com/CGATOxford/UMI-tools), [readthedocs](https://readthedocs.org/projects/umi-tools/) |
| Nextflow | [Nextflow](https://www.nextflow.io/), [Documentation](https://www.nextflow.io/docs/latest/index.html), [GitHub](https://github.com/nextflow-io/nextflow) |

**Note:** Cutadapt v1.18 (2018-09-07), or later, is required.

**Note:** Samtools 1.9 is required if running on CentOS 7.

**Note:** Nextflow 20, or later, is required.

Install:

```console
$ conda install -y pyyaml
$ conda install -y gitpython
$ conda install -y pytest
$ conda install -y pandas
$ conda install -y -c bioconda cutadapt
$ conda install -y -c bioconda pysam
$ conda install -y -c bioconda samtools=1.9
$ conda install -y -c anaconda biopython
$ conda install -y -c bioconda gffutils
$ conda install -y -c anaconda h5py
$ conda install -y -c bioconda umi_tools
$ conda install -y -c bioconda nextflow=20
```

Check packages have installed command-line tools:

```console
$ cutadapt --version
$ samtools --version
$ umi_tools -v
$ javac -version
$ java -version
$ nextflow -v
$ nextflow -version
```

**Note:** Java is installed as a side-effect of installing Nextflow.

Check h5py package has installed:

```console
$ python
```
```python
>>> import h5py
>>> h5py.run_tests()
...
============ 521 passed, 25 skipped, 3 xfailed, 1 warning in 5.50s =============
```

Values for `skipped`, `xfailed` (expected failures) and `warning` may differ, depending upon the version of h5py installed.

Run Nextflow "hello" example [main.nf](https://github.com/nextflow-io/hello/blob/master/main.nf) from [nextflow-io/hello.git](https://github.com/nextflow-io/hello.git):

```console
$ nextflow run hello
N E X T F L O W  ~  version 20.01.0
Pulling nextflow-io/hello ...
downloaded from https://github.com/nextflow-io/hello.git
Launching `nextflow-io/hello` [spontaneous_magritte] - revision: 1d43afc0ec [master]
WARN: The use of `echo` method is deprecated
executor >  local (4)
[1d/bb459e] process > sayHello [100%] 4 of 4 /
Hola world!

Bonjour world!

Ciao world!

Hello world!
```

### Troubleshooting: Samtools `error while loading shared libraries: libcrypto.so.1.0.0`

If you get the following then using Samtools 1.7:

```conda
$ samtools --version
samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory
```

then try forcing a reinstall to a more recent version, for example:

```console
$ conda install -c bioconda samtools=1.9 --force-reinstall
$ samtools --version
samtools 1.9
```

---

## Install tools not available as operating system packages

| Package | Links |
| ------- | ----- |
| Hisat2 (2.1.0) | [Hisat2](https://daehwankimlab.github.io/hisat2/) |
| Bowtie | [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) |

The directory names may differ, depending on the versions you have.

### Install Hisat2

**Note:** Hisat 2.1.0 is strongly recommended. Hisat2 2.2.0 users have reported bugs and issues (see for example [DaehwanKimLab/hisat2#242](https://github.com/DaehwanKimLab/hisat2/issues/242) and [DaehwanKimLab/hisat2#245](https://github.com/DaehwanKimLab/hisat2/issues/245)) which Hisat2 will resolve in a future release.

```console
$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
$ unzip hisat2-2.1.0-Linux_x86_64.zip
$ ls hisat2-2.1.0/
```

Update `PATH` and check that the `hisat2` tool is available:

```console
$ export PATH=~/hisat2-2.1.0:$PATH
$ hisat2 --version
```

### Install Bowtie

```console
$ wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip/download -O bowtie.zip
$ unzip bowtie.zip
$ ls bowtie-1.2.2-linux-x86_64/
```

Update `PATH` and check that the `bowtie` tool is available:

```console
$ export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
$ bowtie --version
```

---

## Create `set-riboviz-env.sh` to configure paths

Create a `set-riboviz-env.sh` script with the paths to your Hisat2 and Bowtie directories. For example:

```console
#!/usr/bin/env bash
export PATH=$HOME/hisat2-2.1.0:$PATH
export PATH=$HOME/bowtie-1.2.2-linux-x86_64/:$PATH
source $HOME/miniconda3/bin/activate
conda activate riboviz
```

Remember, your directory names may differ, depending on the versions of Hisat2 and Bowtie you have.

In future you can configure the paths by running:

```console
$ source set-riboviz-env.sh
```

---

## Get RiboViz

Get RiboViz:

```console
$ git clone https://github.com/riboviz/riboviz
```

---

## Check installation

You can now check your installation by running RiboViz tests by running a "vignette" of the **RiboViz** workflow to see **RiboViz**'s capabilities. See [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md).

Once you have run the "vignette", you can check your installation by running tests:

* [Run vignette regression tests](../developer/testing.md#run-vignette-regression-tests).
* [Run Python tests and workflow tests](../developer/testing.md#run-python-tests-and-workflow-tests).

---

## Reference

### Check names and versions of Python packages

Run:

```console
$ conda list
```

Alternatively, run:

```console
$ pip list
```

The Python packages and their versions will be listed.

### Check names and versions of R packages

Run:

```console
$ Rscript rscripts/list-r-packages.R
```

The R packages and their versions will be listed.
