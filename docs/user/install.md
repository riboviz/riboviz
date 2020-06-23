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
| cURL | 7.68.0 |
| bedtools | 2.26.0 |
| hdf5tools (h5diff) | 1.10.4 |
| pigz | 2.4 |
| Python | 3.7.3 |
| Cutadapt | 1.18 |
| samtools | 1.9 |
| UMI-tools | 1.0.1 |
| Java (javac) | 1.8.0_152-release |
| Java (java) | 1.8.0_152-release |
| Nextflow | 20.01.0.5264 |
| GraphViz (dot) | 2.40.1 |
| Hisat2 | 2.1.0 |
| Bowtie | 1.2.2 |
| R | 3.4.4 |
 
| Python Package | Version | Package Manager |
| -------------- | ------- | --------------- |
| biopython | 1.76 | conda | |
| cutadapt | 1.18 | conda | |
| gitpython | 3.0.5 | conda | |
| h5py | 2.10.0 | conda | |
| nextflow | 20.01.0 | conda | |
| pandas | 1.0.1 | conda | |
| pycodestyle | 2.5.0 | conda | |
| pylint | 2.4.4 | conda | |
| pysam | 0.15.3 | conda | |
| pytest | 5.3.5 | conda | |
| pytest-cov | 2.8.1 | conda | |
| pyyaml | 5.3 | conda | |
| samtools | 1.9 | conda | |
| umi_tools | 1.0.1 | conda | |
| gffutils | 0.10.1 | pip |
| sphinx | 2.4.3 | pip |
 
| R Package | Version |
| --------- | ------- |
| Biostrings | 2.46.0 |
| ggplot2 | 3.2.1 |
| git2r | 0.26.1 |
| here | 0.1 |
| lintr | 2.0.0 |
| optparse | 1.6.4 |
| plotly | 4.9.0 |
| RcppRoll | 0.3.0 |
| readr | 1.3.1 |
| rhdf5 | 2.22.0 |
| Rsamtools | 1.30.0 |
| rtracklayer | 1.38.3 |
| shiny | 1.3.2 |
| tidyr | 1.0.0 |
| ShortRead | 1.36.1 |
| styler | 1.2.0 |
 
Certain packages are only required if you plan to develop and extend RiboViz. These packages are (see [Install developer dependencies](../developer/install.md)):

* Python pycodestyle, pylint, pytest-cov, sphinx
* R: lintr, styler

Constraints:

* RiboViz is **not** compatible with Python 2. Python 2 comes to the end of its supported life in 2020 and there will be no Python 2.8 (see [PEP 373 Python 2.7 Release Schedule](https://legacy.python.org/dev/peps/pep-0373/)).
* Either [Miniconda](https://conda.io/miniconda.html) Python 3.6, or later, or [Anaconda Distribution](https://www.anaconda.com/distribution/) Python 3.6, or later, are strongly recommended.
* Cutadapt v1.18 (2018-09-07), or later, is required.
* R 2.14.0, or later, is required as it includes the [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html) package.
* R 3.4, or later, is strongly recommended.

---

## General packages

### Git

Website: [Git](https://git-scm.com/)

**Ubuntu**

Install:

```console
$ sudo apt-get install -y git
```

**CentOS**

Install:

```console
$ sudo yum install -y git
```

**Check package has installed:**

```console
$ git --version
```

### cURL

Website: [cURL](https://curl.haxx.se/)

**Ubuntu**

Install:

```console
$ sudo apt-get install -y curl
```

**CentOS**

Install:

```console
$ sudo yum install -y curl
```

**Check package has installed:**

```console
$ curl --version
```

### EPEL (CentOS only)

Website [EPEL](https://fedoraproject.org/wiki/EPEL)

Install:

```console
$ sudo yum install -y epel-release
```

### bedtools

Web sites:

* [bedtools](http://bedtools.readthedocs.io/en/latest/)
* [GitHub](https://github.com/arq5x/bedtools2)

**Ubuntu**

Install:

```console
$ sudo apt-get install -y bedtools
```

**CentOS**

Install:

```console
$ sudo yum install -y BEDTools
```

**Check package has installed:**

```console
$ bedtools -version
```

### hdf5tools

Web site: [HDF5](https://portal.hdfgroup.org/display/HDF5)

**Ubuntu**

Install:

```console
$ sudo apt-get install -y hdf5-tools
```

**CentOS**

Install:

```console
$ sudo yum install -y hdf5-devel
```

**Check package has installed:**

```console
$ h5diff -version
```

`h5diff` is one of the hdf5tools.

### pigz

Web site: [pigz](http://zlib.net/pigz/)

pigz is a parallel implementation of gzip for multi-processor, multi-core machines. It may already be available on your system. Check by running the following command:

```console
$ pigz --version
```

If it is not present, please install it as follows:

**Ubuntu**

Install:

```console
$ sudo apt-get install -y pigz
```

**CentOS**

Install:

```console
$ sudo yum install -y pigz
```

---

## Python

Web site: [python](https://www.python.org/)

RiboViz is **not** compatible with Python 2. Python 2 comes to the end of its supported life in 2020 and there will be no Python 2.8 (see [PEP 373 Python 2.7 Release Schedule](https://legacy.python.org/dev/peps/pep-0373/)).

Either [Miniconda](https://conda.io/miniconda.html) Python 3.6, or later, or [Anaconda Distribution](https://www.anaconda.com/distribution/) Python 3.6, or later, are strongly recommended.

The instructions which follow have been written under the assumption that you are using Miniconda Python. If using Anaconda then, when installing some packages, you will be told that they are already available. This is because Anaconda comes with a wide range of common Python packages.

If you are using other distributions of Python you will need to consult the relevant documentation for each package for installation information.

### Miniconda Python 3.6+

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

## Python packages

### pyyaml

Web sites:

* [PyYAML](https://pyyaml.org/)
* [GitHub](https://github.com/yaml/pyyaml/)

Install:

```console
$ conda install -y pyyaml
```

### gitpython

Web sites:

* [gitpython](https://gitpython.readthedocs.io/en/stable/)
* [GitHub](https://github.com/gitpython-developers/GitPython)

Install:

```console
$ conda install -y gitpython
```

### pytest

Web sites:

* [pytest](https://pytest.org/)
* [GitHub](https://github.com/pytest-dev/pytest/)

Install:

```console
$ conda install -y pytest
```

### pandas

Web sites:

* [pandas](https://pandas.pydata.org/)
* [GitHub](https://github.com/pandas-dev/pandas)

Install:

```console
$ conda install -y pandas
```

### Cutadapt

Web sites:

* [GitHub](https://github.com/marcelm/cutadapt)
* [readthedocs](https://cutadapt.readthedocs.io/)

Cutadapt v1.18 (2018-09-07), or later, is required.

Install:

```console
$ conda install -y -c bioconda cutadapt
```

Check package has installed the `cutadapt` tool:

```console
$ cutadapt --version
```

### pysam

Web sites:

* [GitHub](https://github.com/pysam-developers/pysam/)
* [readthedocs](https://pysam.readthedocs.io/)

Install:

```console
$ conda install -y -c bioconda pysam
$ conda install -y -c bioconda samtools
```

Check package has installed the `samtools` tool:

```console
$ samtools --version
```

### BioPython

Web site:

* [Biopython](http://biopython.org/)

Install:

```console
$ conda install -y -c anaconda biopython
```

### gffutils

Web site:

* [gffutils](http://daler.github.io/gffutils/)

Install:

```console
$ pip install gffutils
```

**Note:** `pip install` is recommended. Using

```console
$ conda install -y -c bioconda gffutils
```

under Python 3, seems to confuse the Python environment and sets Python to:

```console
$ python --version
Python 2.7.16 :: Anaconda, Inc.
```

### h5py

Web site:

* [h5py](https://www.h5py.org/)

Install:

```console
$ conda install -y -c anaconda h5py
```

Check package has installed:

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

### UMI-tools

Web sites:

* [GitHub](https://github.com/CGATOxford/UMI-tools)
* [readthedocs](https://readthedocs.org/projects/umi-tools/)

Install:

```console
$ conda install -y -c bioconda umi_tools
```

Check package has installed the `umi_tools` tool:

```console
$ umi_tools -v
```

---

## Hisat2

Web site: [Hisat2](https://daehwankimlab.github.io/hisat2/)

You may choose to install a more recent version of Hisat2 than in the example that follows.

Install:

```console
$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
$ unzip hisat2-2.1.0-Linux_x86_64.zip
$ ls hisat2-2.1.0/
```

Your directory name may differ, depending on the version of Hisat2 you have.

Update `PATH` and check that the `hisat2` tool is available:

```console
$ export PATH=~/hisat2-2.1.0:$PATH
$ hisat2 --version
```

---

## Bowtie

Web site: [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)

**Note:** We are working to add back a Bowtie alignment option.

You may choose to install a more recent version of Hisat2 than in the example that follows.

Install:

```console
$ wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip/download -O bowtie.zip
$ unzip bowtie.zip
$ ls bowtie-1.2.2-linux-x86_64/
```

Your directory name may differ, depending on the version of Bowtie you have.

Update `PATH` and check that the `bowtie` tool is available:

```console
$ export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
$ bowtie --version
```

## Create `setenv.sh` to configure Hisat2 and Bowtie paths

Create a `setenv.sh` script with the paths to your Hisat2 and Bowtie directories. For example:

```console
#!/usr/bin/env bash
export PATH=~/hisat2-2.1.0:$PATH
export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
```

Remember, your directory names may differ, depending on the versions of Hisat2 and Bowtie you have.

In future you can configure the paths by running:

```console
$ source setenv.sh
```

---

## R 2.14.0+

Web sites:

* [The R Project for Statistical Computing](https://www.r-project.org/)
* [The Comprehensive R Archive Network](https://cran.r-project.org/) (CRAN).

R 2.14.0, or later, is required as it includes the [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html) package.

R 3.4, or later, is strongly recommended.

### Install R and packages required by R packages to be installed

**Ubuntu**

```console
$ sudo apt-get update -y
$ sudo apt-get install -y r-base
$ sudo apt-get install -y r-base-dev

$ sudo apt-get install -y libxml2-dev
$ sudo apt-get install -y libssl-dev
$ sudo apt-get install -y libcurl4-openssl-dev
```

**CentOS**

```console
$ sudo yum update -y
$ sudo yum install -y R
$ sudo yum install -y R-devel

$ sudo yum install -y libxml2-devel
$ sudo yum install -y openssl-devel
$ sudo yum install -y libcurl-devel
```

**Troubleshooting: the most recent version of R is not installed**

Default package managers may not have the most up-to-date version of R available. [The Comprehensive R Archive Network](https://cran.r-project.org/) has information on alternative ways to get a more recent version of R.

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
$ sudo apt-get -s install r-base
$ R --version
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
```

### Check R has installed

```console
$ R --version
R version 3.4.4 (2018-03-15) -- "Someone to Lean On"
$ Rscript --version
R scripting front-end version 3.4.4 (2018-03-15)
```

Your version of R may differ from that shown.

---

## R packages

If, when installing R packages, you see a message like:

```
Would you like to use a personal library instead? (yes/No/cancel) yes
Would you like to create a personal library
"~/R/x86_64-pc-linux-gnu-library/3.6"
to install packages into? (yes/No/cancel)
```

then enter `yes`.

### RcppRoll

Web site: [RcppRoll](https://cran.r-project.org/web/packages/RcppRoll/index.html)

Install in R:

```r
> install.packages("RcppRoll")
```

### optparse

Web site: [optparse](https://cran.r-project.org/web/packages/optparse/index.html)

Install in R:

```r
> install.packages("optparse")
```

### tidyr

Web site: [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)

Install in R:

```r
> install.packages("tidyr")
```

### ggplot2

Web site: [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)

Install in R:

```r
> install.packages("ggplot2")
```

### shiny

Web site: [shiny](https://cran.r-project.org/web/packages/shiny/index.html)

Install in R:

```r
> install.packages("shiny")
```

### plotly

Web site: [plotly](https://cran.r-project.org/web/packages/plotly/index.html)

Install in R:

```r
> install.packages("plotly")
```

### readr

Website: [readr](https://cran.r-project.org/web/packages/readr/index.html)

Install in R:

```r
> install.packages("readr")
```

### git2r

Web sites:

* [git2r](https://docs.ropensci.org/git2r)
* [GitHub](https://github.com/ropensci/git2r)

Install in R:

```r
> install.packages("git2r")
```

### here

Web sites:

* [here](https://here.r-lib.org/)
* [CRAN](https://cran.r-project.org/package=here)
* [GitHub](https://github.com/r-lib/here)

Install in R:

```r
> install.packages("here")
```

### Bioconductor Rsamtools, rhdf5, rtracklayer, Biostrings, ShortRead

Web sites:

* [Bioconductor](https://bioconductor.org)
* [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html)
* [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
* [rhdf5](https://bioconductor.org/packages/release/bioc/html/rhdf5.html)
* [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
* [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html)

The commands to install Bioconductor packages depend on your version of R. For full details:

* See Bioconductor/R compatibility on [Bioconductor releases](https://bioconductor.org/about/release-announcements/).
* Click the link of a Bioconductor release consistent with your version of R.
* Click the link of the specific package.

For example, for R 3.4, install in R:

```r
> source("https://bioconductor.org/biocLite.R")
> biocLite("Rsamtools")
> biocLite("rtracklayer")
> biocLite("rhdf5")
> biocLite("Biostrings")
> biocLite("ShortRead")
```

For example, for R 3.5, install in R:

```r
> install.packages("BiocManager")
> BiocManager::install("Rsamtools")
> BiocManager::install("rtracklayer")
> BiocManager::install("rhdf5")
> BiocManager::install("Biostrings")
> BiocManager::install("ShortRead")
```

**Troubleshooting: installation path not writeable**

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

**Troubleshooting: Cannot allocate memory**

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

---

## Nextflow

**Note:** Nextflow only needs to be installed if you plan to use the Nextflow version of RiboViz.

Web sites:

* [Nextflow](https://www.nextflow.io/)
* [Documentation](https://www.nextflow.io/docs/latest/index.html)
* [GitHub](https://github.com/nextflow-io/nextflow)

### Install Nextflow using conda (recommended)

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

### Install Nextflow (alternative)

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

Update the `setenv.sh` script (see [Create `setenv.sh` to configure Hisat2 and Bowtie paths](#create-setenvsh-to-configure-hisat2-and-bowtie-paths) above) by adding:

```
export PATH=$HOME/nextflow:$PATH
```

### Run Nextflow "hello" example

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

## RiboViz

Get RiboViz:

```console
$ git clone https://github.com/riboviz/riboviz
```

---

## Check installation

You can now check your installation by running RiboViz tests.

Run tests:

```console
$ pytest --ignore-glob="*regression*" --ignore-glob="*nextflow*"
=============================== test session starts ===============================
platform linux -- Python 3.7.3, pytest-5.3.5, py-1.8.1, pluggy-0.13.1
rootdir: /home/ubuntu/riboviz
plugins: cov-2.8.1
collected 166 items                                                               
riboviz/test/test_barcodes_umis.py .................                        [ 10%]
riboviz/test/test_demultiplex_fastq.py ..................                   [ 21%]
riboviz/test/test_fastq.py ........................                         [ 35%]
riboviz/test/test_process_utils.py ........................                 [ 50%]
riboviz/test/test_sam_bam.py ................                               [ 59%]
riboviz/test/test_trim_5p_mismatch.py ............                          [ 66%]
riboviz/test/test_upgrade_config.py .....                                   [ 69%]
riboviz/test/test_utils.py ....................                             [ 81%]
riboviz/test/tools/test_prep_riboviz.py ............                        [ 89%]
riboviz/test/tools/test_prep_riboviz_simdata_multiplex.py 
..............    [ 97%]
riboviz/test/tools/test_prep_riboviz_simdata_umi.py ....                    [100%]

...
=================== 166 passed, 1 warning in 259.81s (0:04:19) ====================
```

`PendingDeprecationWarning` `warnings` can be ignored.

If you installed Nextflow, run the Nextflow tests too:

```console
$ pytest riboviz/test/nextflow
============================== test session starts ===============================
platform linux -- Python 3.7.3, pytest-5.3.5, py-1.8.1, pluggy-0.13.1
rootdir: /home/ubuntu/riboviz
plugins: cov-2.8.1
collected 36 items

riboviz/test/nextflow/test_nextflow_errors.py ............................ [ 77%]
........                                                                   [100%]

========================= 36 passed in 143.72s (0:02:23) =========================
```

Download regression test data:

```console
$ cd
$ git clone https://github.com/riboviz/regression-test-data-2.0.beta
```

Run the regression tests for the RiboViz Python workflow (these may take a few minutes):

```console
$ pytest riboviz/test/regression/test_regression.py --expected=$HOME/regression-test-data-2.0.beta/
============================= test session starts ==============================
platform linux -- Python 3.7.3, pytest-5.3.5, py-1.8.1, pluggy-0.13.1 -- /home/ubuntu/miniconda3/bin/python
rootdir: /home/ubuntu/riboviz
collected 82 items

riboviz/test/regression/test_regression.py ssssssssssssssssssssssssssssssss [ 39%]
ssssss..ssss......................................                          [100%]

...
============== 40 passed, 42 skipped, 1 warning in 215.47s (0:03:35) ==============
```

If you installed Nextflow, run the regression tests for the RiboViz Nextflow workflow (these may take a few minutes):

```console
$ pytest riboviz/test/regression/test_regression.py --expected=$HOME/regression-test-data-2.0.beta/ --nextflow
=============================== test session starts ===============================
platform linux -- Python 3.7.3, pytest-5.3.5, py-1.8.1, pluggy-0.13.1 -- /home/ubuntu/miniconda3/bin/python
cachedir: .pytest_cache
rootdir: /home/ubuntu/riboviz
plugins: cov-2.8.1
collected 82 items

riboviz/test/regression/test_regression.py ssssssssssssssssssssssssssssssss [ 39%]
ssssss..ssss......................................                          [100%]

============== 40 passed, 42 skipped, 1 warning in 163.16s (0:02:43) ==============
```

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
