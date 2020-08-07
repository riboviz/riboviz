# Install RiboViz and dependencies

## About these instructions

These instructions were written for Ubuntu 18.04 and CentOS 7.4. Other Linux flavours will require different commands to be run.

Installing some of these tools requires you to have permission to run `sudo` to install and configure software. If you don't have `sudo` access you will have to ask a local system administrator to run these commands for you.

### Windows users

We suggest that you:

* Either, use a virtual machine running under [VMWare Workstation Player](https://www.vmware.com/uk/products/workstation-player.html) or [Oracle VirtualBox](https://www.virtualbox.org/).
* Or, try using Windows 10's [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) which allows running of a Linux environment on Windows 10 without the need for a VM. Most command-line tools, utilities, and applications can be directly on Windows, unmodified. Ubuntu, openSUSE, Debian, Kali flavours of Linux can be used.

### Mac OSX users

We suggest that you check out the web sites for each prerequisite for information on how to install the prerequisites under Mac OS X.

---

## Dependencies overview

The following tables summarise the packages required by RiboViz. Instructions to install each dependency are given in the following sections - only minimal installation instructions are given, see the documentation for each dependency for full instructions.

The versions listed are those used by a RiboViz developer when preparing the current release. Other versions of the prerequisites, different from the versions shown, may also be usable, but see the constraints below.

| Command-line tool | Version |
| ----------------- | ------- |
| Git | 2.17.1 |
| cURL | 7.71.0 |
| bedtools | 2.26.0 |
| hdf5tools (h5diff) | 1.10.6 |
| pigz | 2.4 |
| Python | 3.7.7 |
| Cutadapt | 1.18 |
| samtools | 1.7 |
| UMI-tools | 1.0.1 |
| Java (javac) | 1.8.0_152-release |
| Java (java) | 1.8.0_152-release |
| Nextflow | 20.04.1.5335 |
| GraphViz (dot) | 2.40.1 |
| Hisat2 | 2.1.0 |
| Bowtie | 1.2.2 |
| R | 3.6.3 |
 
| Python Package | Version | Package Manager |
| -------------- | ------- | --------------- |
| biopython | 1.77 | conda | |
| cutadapt | 1.18 | conda | |
| gitpython | 3.1.3 | conda | |
| h5py | 2.10.0 | conda | |
| nextflow | 20.04.1 | conda | |
| pandas | 1.0.5 | conda | |
| pycodestyle | 2.6.0 | conda | |
| pylint | 2.5.3 | conda | |
| pysam | 0.15.3 | conda | |
| pytest | 5.4.3 | conda | |
| pytest-cov |  | conda | |
| pyyaml | 5.3.1 | conda | |
| samtools | 1.7 | conda | |
| umi_tools | 1.0.1 | conda | |
| gffutils | 0.10.1 | pip |
| sphinx |  | pip |
 
| R Package | Version |
| --------- | ------- |
| Biostrings | 2.54.0 |
| ggplot2 | 3.3.2 |
| git2r | 0.27.1 |
| here | 0.1 |
| lintr | 2.0.1 |
| optparse | 1.6.6 |
| plotly | 4.9.2.1 |
| RcppRoll | 0.3.0 |
| readr | 1.3.1 |
| rhdf5 | 2.30.1 |
| Rsamtools | 2.2.3 |
| rtracklayer | 1.46.0 |
| shiny | 1.5.0 |
| tidyr | 1.1.0 |
| ShortRead | 1.44.3 |
| styler | 1.3.2 |

Certain packages are only required if you plan to develop and extend RiboViz. These packages are (see [Install developer dependencies](../developer/install.md)):

* Python pycodestyle, pylint, pytest-cov, sphinx
* R: lintr, styler

Constraints:

* RiboViz is **not** compatible with Python 2. Python 2 comes to the end of its supported life in 2020 and there will be no Python 2.8 (see [PEP 373 Python 2.7 Release Schedule](https://legacy.python.org/dev/peps/pep-0373/)).
* Either [Miniconda](https://conda.io/miniconda.html) Python 3.6, or later, or [Anaconda Distribution](https://www.anaconda.com/distribution/) Python 3.6, or later, are strongly recommended.
* Cutadapt v1.18 (2018-09-07), or later, is required.
* Hisat 2.1.0 is recommended, not 2.2.0. Hisat2 2.2.0 users have reported bugs and issues (see for example [DaehwanKimLab/hisat2#242](https://github.com/DaehwanKimLab/hisat2/issues/242) and [DaehwanKimLab/hisat2#245](https://github.com/DaehwanKimLab/hisat2/issues/245)) which Hisat2 say will be resolved in a future release.
* R 2.14.0, or later, is required as it includes the [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html) package.
* R 3.6, or later, is strongly recommended.

---

## Install general packages

| Package | Links |
| ------- | ----- |
| Git | [Git](https://git-scm.com/) |
| cURL | [cURL](https://curl.haxx.se/) |
| EPEL (CentOS only) | [EPEL](https://fedoraproject.org/wiki/EPEL) |
| bedtools | [bedtools](http://bedtools.readthedocs.io/en/latest/), [GitHub](https://github.com/arq5x/bedtools2) |
| hdf5tools | [HDF5](https://portal.hdfgroup.org/display/HDF5) |
| pigz | [pigz](http://zlib.net/pigz/) |

### Install on Ubuntu

```console
$ sudo apt-get install -y git
$ sudo apt-get install -y curl
$ sudo apt-get install -y bedtools
$ sudo apt-get install -y hdf5-tools
$ sudo apt-get install -y pigz
```

### Install on CentOS

```console
$ sudo yum install -y git
$ sudo yum install -y curl
$ sudo yum install -y epel-release
$ sudo yum install -y BEDTools
$ sudo yum install -y hdf5-devel
$ sudo yum install -y pigz
```

### Check packages have installed

```console
$ git --version
$ curl --version
$ bedtools -version
$ h5diff -version
$ pigz --version
```

`h5diff` is one of the hdf5tools.

---

## Install Python

Web site: [python](https://www.python.org/)

RiboViz is **not** compatible with Python 2. Python 2 comes to the end of its supported life in 2020 and there will be no Python 2.8 (see [PEP 373 Python 2.7 Release Schedule](https://legacy.python.org/dev/peps/pep-0373/)).

Either [Miniconda](https://conda.io/miniconda.html) Python 3.6, or later, or [Anaconda Distribution](https://www.anaconda.com/distribution/) Python 3.6, or later, are strongly recommended.

The instructions which follow have been written under the assumption that you are using Miniconda Python. If using Anaconda then, when installing some packages, you will be told that they are already available. This is because Anaconda comes with a wide range of common Python packages.

If you are using other distributions of Python you will need to consult the relevant documentation for each package for installation information.

### Install Miniconda Python 3.6+

Install:

```console
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
$ bash miniconda3.sh -b -p $HOME/miniconda3
```

**Note:** make sure you use `-O`, which provides a name for the downloaded file, and not `-o`, which provides the name of a file for messages about the download.

Activate environment:

```console
$ source $HOME/miniconda3/bin/activate
```

Create a `riboviz` environment and activate it:

```console
$ conda create --name riboviz python=3.7
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

**Troubleshooting: incompatible Python versions**

If you find when installing packages below that your version of Python is too new (e.g. your Python version is 3.8 and the tool only works with Python 3.7), then two options are:

1. See if a more recent version of the tool is available. For example a version available via `pip` may be more up-to-date than a version available via `conda`.
2. Create a conda environment that supports a version of Python compatible with the tools.

```console
$ conda create --name riboviz python=<VERSION>
```

  - For example:

```console
$ conda create --name riboviz python=3.7
```

---

## Install Python packages

| Package | Links |
| ------- | ----- |
| pyyaml | [PyYAML](https://pyyaml.org/), [GitHub](https://github.com/yaml/pyyaml/) |
| gitpython | [gitpython](https://gitpython.readthedocs.io/en/stable/), [GitHub](https://github.com/gitpython-developers/GitPython) |
| pytest | [pytest](https://pytest.org/), [GitHub](https://github.com/pytest-dev/pytest/) |
| pandas | [pandas](https://pandas.pydata.org/), [GitHub](https://github.com/pandas-dev/pandas) |
| Cutadapt | [GitHub](https://github.com/marcelm/cutadapt), [readthedocs](https://cutadapt.readthedocs.io/) |
| pysam | [GitHub](https://github.com/pysam-developers/pysam/), [readthedocs](https://pysam.readthedocs.io/) |
| BioPython | [Biopython](http://biopython.org/) |
| gffutils | [gffutils](http://daler.github.io/gffutils/) |
| h5py | [h5py](https://www.h5py.org/) |
| UMI-tools | [GitHub](https://github.com/CGATOxford/UMI-tools), [readthedocs](https://readthedocs.org/projects/umi-tools/) |

**Note:** Cutadapt v1.18 (2018-09-07), or later, is required.

Install:

```console
$ conda install -y pyyaml
$ conda install -y gitpython
$ conda install -y pytest
$ conda install -y pandas
$ conda install -y -c bioconda cutadapt
$ conda install -y -c bioconda pysam
$ conda install -y -c bioconda samtools
$ conda install -y -c anaconda biopython
$ pip install gffutils
$ conda install -y -c anaconda h5py
$ conda install -y -c bioconda umi_tools
```

Check packages have installed command-line tools:

```console
$ cutadapt --version
$ samtools --version
$ umi_tools -v
```

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

Your number of `skipped` and `xfailed` (expected failures may differ, depending upon the version of h5py installed.

**Note:** For `gffutils`, `pip install` is recommended because using

```console
$ conda install -y -c bioconda gffutils
```

under Python 3, seems to confuse the Python environment and sets Python to:

```console
$ python --version
Python 2.7.16 :: Anaconda, Inc.
```

---

## Install Bioinformatics tools

| Package | Links |
| ------- | ----- |
| Hisat2 (2.1.0) | [Hisat2](https://daehwankimlab.github.io/hisat2/) |
| Bowtie | [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) |
| Nextflow | [Nextflow](https://www.nextflow.io/), [Documentation](https://www.nextflow.io/docs/latest/index.html), [GitHub](https://github.com/nextflow-io/nextflow) |

The directory names may differ, depending on the versions you have.

### Install Hisat2

**Note:** Hisat 2.1.0 is recommended, not 2.2.0. Hisat2 2.2.0 users have reported bugs and issues (see for example [DaehwanKimLab/hisat2#242](https://github.com/DaehwanKimLab/hisat2/issues/242) and [DaehwanKimLab/hisat2#245](https://github.com/DaehwanKimLab/hisat2/issues/245)) which Hisat2 say will be resolved in a future release.

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

**Note:** We are working to add a Bowtie alignment option.

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

### Install Nextflow

**Note:** Nextflow only needs to be installed if you plan to use the Nextflow version of RiboViz.

**Install Nextflow using conda (recommended)**

Install Nextflow and its dependencies (including Java):

```console
$ conda install -y -c bioconda nextflow
```

Check install:

```console
$ javac -version
$ java -version
$ nextflow -version

      N E X T F L O W
      version 20.01.0 build 5264
      created 12-02-2020 10:14 UTC (02:14 PDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```


Your version of Nextflow may differ from that shown.

**Install Nextflow (alternative)**

Install [OpenJDK](https://openjdk.java.net) 1.8:

* Ubuntu 18 users:

```console
$ sudo apt-get install -y openjdk-8-jdk-headless
```

* CentOS 7 users:

```console
$ sudo yum install -y openjdk-8-jdk-headless
```

Check install:

```console
$ javac -version
$ java -version
```

Install Nextflow:

```console
$ curl -s https://get.nextflow.io | bash
$ export PATH=$HOME/nextflow:$PATH
$ nextflow -version
```

Set `PATH`:

```console
$ export PATH=$HOME/nextflow:$PATH
```

**Run Nextflow "hello" example**

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

This runs Nextflow workflow [main.nf](https://github.com/nextflow-io/hello/blob/master/main.nf) from [nextflow-io/hello.git](https://github.com/nextflow-io/hello.git).

---

## Create `setenv.sh` to configure paths

Create a `setenv.sh` script with the paths to your Hisat2 and Bowtie directories. For example:

```console
#!/usr/bin/env bash
export PATH=~/hisat2-2.1.0:$PATH
export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
```

Remember, your directory names may differ, depending on the versions of Hisat2 and Bowtie you have.

If you installed Nextflow using the "Install Nextflow (alternative)" instructions then also add:

```
export PATH=$HOME/nextflow:$PATH
```

In future you can configure the paths by running:

```console
$ source setenv.sh
```

---

## Install R 2.14.0+

Web sites:

* [The R Project for Statistical Computing](https://www.r-project.org/)
* [The Comprehensive R Archive Network](https://cran.r-project.org/) (CRAN).

**Note:** R 2.14.0, or later, is required as it includes the [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html) package.

**Note:** R 3.6, or later, is strongly recommended.

### Install R and packages required by R packages to be installed

**Install on Ubuntu**

```console
$ sudo apt-get update -y
$ sudo apt-get install -y r-base
$ sudo apt-get install -y r-base-dev

$ sudo apt-get install -y libxml2-dev
$ sudo apt-get install -y libssl-dev
$ sudo apt-get install -y libcurl4-openssl-dev
```

**Install on CentOS**

```console
$ sudo yum update -y
$ sudo yum install -y R
$ sudo yum install -y R-devel

$ sudo yum install -y libxml2-devel
$ sudo yum install -y openssl-devel
$ sudo yum install -y libcurl-devel
```

**Troubleshooting: `the most recent version of R is not installed` or `package "..." is not available (for R version ...)`

Default package managers may not have the most up-to-date version of R available. This may cause problems when installing R packages. [The Comprehensive R Archive Network](https://cran.r-project.org/) has information on alternative ways to get a more recent version of R.

For example, following CRAN-R's [UBUNTU PACKAGES FOR R](https://cran.r-project.org/bin/linux/ubuntu/README.html) to get the latest R 3.6 packages:

* Get your Ubuntu version, for example:

```console
$ lsb_release -a
No LSB modules are available.
Distributor ID:	Ubuntu
Description:	Ubuntu 18.04 LTS
Release:	18.04
Codename:	bionic
```

* Open `/etc/apt/sources.list` in an editor (requires `sudo` access), for example:

```console
$ sudo nano /etc/apt/sources.list
```

* From UBUNTU PACKAGES FOR R get the entry for your Ubuntu version and add it to the end of the file. For example:

```
deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/
```

* Install CRAN-R Ubuntu server key:

```console
$ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
```

* Update packages:

```console
$ sudo apt-get update
```

* Install:

```console
$ sudo apt-get install -y r-base
$ sudo apt-get install -y r-base-dev
$ R --version
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
```

### Check R has installed

```console
$ R --version
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
```

Your version of R may differ from that shown.

---

## Install R packages

| Package | Links |
| ------- | ----- |
| RcppRoll | [RcppRoll](https://cran.r-project.org/web/packages/RcppRoll/index.html) |
| optparse | [optparse](https://cran.r-project.org/web/packages/optparse/index.html) |
| tidyr | [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html) |
| ggplot2 | [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) |
| shiny | [shiny](https://cran.r-project.org/web/packages/shiny/index.html) |
| plotly | [plotly](https://cran.r-project.org/web/packages/plotly/index.html) |
| readr | [readr](https://cran.r-project.org/web/packages/readr/index.html) |
| git2r |  [git2r](https://docs.ropensci.org/git2r), [GitHub](https://github.com/ropensci/git2r) |
| here | [here](https://here.r-lib.org/), [CRAN](https://cran.r-project.org/package=here), [GitHub](https://github.com/r-lib/here) |
| Bioconductor Rsamtools | [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) |
| Bioconductor rhdf5 | [rhdf5](https://bioconductor.org/packages/release/bioc/html/rhdf5.html) |
| Bioconductor rtracklayer | [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) |
| Bioconductor Biostrings | [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) |
| Bioconductor ShortRead | [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) |

If, when installing R packages, you see a message like:

```
Would you like to use a personal library instead? (yes/No/cancel) yes
Would you like to create a personal library
"~/R/x86_64-pc-linux-gnu-library/3.6"
to install packages into? (yes/No/cancel)
```

then enter `yes`.

Install in R:

```r
> install.packages("RcppRoll")
> install.packages("optparse")
> install.packages("tidyr")
> install.packages("ggplot2")
> install.packages("shiny")
> install.packages("plotly")
> install.packages("readr")
> install.packages("git2r")
> install.packages("here")
```

The commands to install Bioconductor packages depend on your version of R. For full details:

* See Bioconductor/R compatibility on [Bioconductor releases](https://bioconductor.org/about/release-announcements/).
* Click the link of a Bioconductor release consistent with your version of R.
* Click the link of the specific package.

For example, for R 3.5 or R 3.6, install in R:

```r
> install.packages("BiocManager")
> BiocManager::install("Rsamtools")
> BiocManager::install("rtracklayer")
> BiocManager::install("rhdf5")
> BiocManager::install("Biostrings")
> BiocManager::install("ShortRead")
```

For example, for R 3.4, install in R:

```r
> source("https://bioconductor.org/biocLite.R")
> biocLite("Rsamtools")
> biocLite("rtracklayer")
> biocLite("rhdf5")
> biocLite("Biostrings")
> biocLite("ShortRead")
```

### Troubleshooting: installation path not writeable

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
$ ls ~/R/x86_64-redhat-linux-gnu-library/3.5/Rsamtools/
DESCRIPTION  libs  LICENSE  Meta  NAMESPACE  NEWS
```

### Troubleshooting: Cannot allocate memory

```
Error in system2(file.path(R.home("bin"), "R"), c(if (nzchar(arch))
paste0("--arch=",  :
  cannot popen ' '/usr/lib/R/bin/R' --no-save --slave 2>&1 <
  '/tmp/Rtmpw3pOH7/file12471113d0d2b'', probable reason 'Cannot
  allocate memory'
* removing "/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.5/Rsamtools"
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

### Troubleshooting: package "XML" is not available (for R version 3.6.3) 

If you get this errort when running:

```R
> BiocManager::install("rtracklayer")
```

then one solution may be to install "XML" specifying the URL of the source package, for example:

```R
> install.packages("https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz", repos=NULL, type="source")
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

You can now check your installation by running RiboViz tests.

Run tests:

```console
$ cd riboviz
$ pytest --ignore-glob="*regression*" --ignore-glob="*nextflow*"
```

All tests should pass (some may be skipped, but none should fail). `PendingDeprecationWarning` `warnings` can be ignored.

If you installed Nextflow, run the Nextflow tests too:

```console
$ pytest riboviz/test/nextflow
```

Again, all tests should pass (some may be skipped, but none should fail).

Download regression test data:

```console
$ cd
$ git clone https://github.com/riboviz/regression-test-data-2.0
```

Run the regression tests for the RiboViz Python workflow (these may take a few minutes):

```console
$ cd riboviz
$ pytest riboviz/test/regression/test_regression.py --expected=$HOME/regression-test-data-2.0/
```

All tests should pass (some may be skipped, but none should fail).

If you installed Nextflow, run the regression tests for the RiboViz Nextflow workflow (these may take a few minutes):

```console
$ pytest riboviz/test/regression/test_regression.py --expected=$HOME/regression-test-data-2.0/ --nextflow
```

Again, all tests should pass (some may be skipped, but none should fail).

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

`nextflow` will only be shown if you installed Nextflow.

### Check names and versions of R packages

Run:

```console
$ Rscript rscripts/list-r-packages.R
```

The R packages and their versions will be listed.
