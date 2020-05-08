# Install prerequisites

## About these instructions

These instructions were written for Ubuntu 18.04 and CentOS 7.4. Other Linux flavours will require different commands to be run.

Other versions of the prerequisites, different from the versions stated, may be usable

Only minimal installation instructions are given for each prerequisite. See the documentation for each prerequisite for full instructions.

Installing some of these tools requires you to have permission to run `sudo` to install and configure software. If you don't have `sudo` access you will have to ask a local system administrator to run these commands for you.

### Windows users

We suggest that you:

* Either, use a virtual machine running under [VMWare Workstation Player](https://www.vmware.com/uk/products/workstation-player.html) or [Oracle VirtualBox](https://www.virtualbox.org/).
* Or, try using Windows 10's [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) which allows running of a Linux environment on Windows 10 without the need for a VM. Most command-line tools, utilities, and applications can be directly on Windows, unmodified. Ubuntu, openSUSE, Debian, Kali flavours of Linux can be used.

### Mac OSX users

We suggest that you check out the web sites for each prerequisite for information on how to install the prerequisites under Mac OS X.

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
bedtools v2.26.0
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

### pigz

Web site: [pigz](http://zlib.net/pigz/)

pigz is a parallel implementation of gzip for multi-processor, multi-core machines. It may already be available on your system. Check by running the following command:

```console
$ pigz --version
pigz 2.4
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

It is strongly recommended that you use [Miniconda](https://conda.io/miniconda.html) Python 3.6+. Alternatively, use the [Anaconda Distribution](https://www.anaconda.com/distribution/) of Python 3.6+.

The instructions which follow have been written under the assumption that you are using Miniconda Python. If using Anaconda then, when installing some packages, you will be told that they are already available. This is because Anaconda comes with a wide range of common Python packages.

If you are using other distributions of Python you will need to consult the relevant documentation for each package for installation information.

**Note:** RiboViz is **not** compatible with Python 2. Python 2 comes to the end of its supported life in 2020 and there will be no Python 2.8 (see [PEP 373 Python 2.7 Release Schedule](https://legacy.python.org/dev/peps/pep-0373/)).

### Miniconda Python 3.6+

Install:

```console
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
$ bash miniconda3.sh -b -p $HOME/miniconda3
```

**Note:** make sure you use `-O`, which provides a name for the downloaded file, and not `-o`, which provides the name of a file for messages about the download.

Activate environment and check:

```console
$ source $HOME/miniconda3/bin/activate
$ python -V
Python 3.7.3
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

Install:

```console
$ conda install -y -c bioconda cutadapt
```

Check package has installed `cutadapt` tool:

```console
$ cutadapt --v
2.3
```

### pysam

Web sites:

* [GitHub](https://github.com/pysam-developers/pysam/)
* [readthedocs](https://pysam.readthedocs.io/)

Install:

```console
$ conda install -y -c bioconda pysam
```

Check package has installed `samtools` tool:

```console
$ samtools --version
samtools 1.9
Using htslib 1.9
Copyright (C) 2018 Genome Research Ltd.
```

If `samtools` was not installed, then install it explicitly:

```console
$ conda install -c bioconda samtools
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
OK (skipped=16, expected failures=6)
```

The test report should be `OK` as above. Your number of `skipped` and `expected failures` may differ, depending upon the version of h5py installed.

### UMI-tools

Web sites:

* [GitHub](https://github.com/CGATOxford/UMI-tools)
* [readthedocs](https://readthedocs.org/projects/umi-tools/)

Install:

```console
$ conda install -y -c bioconda umi_tools
```

Check package has installed `umi_tools` tool:

```console
$ umi_tools -v
UMI-tools version: 1.0.0
```

---

## Hisat2

Web site: [Hisat2](https://daehwankimlab.github.io/hisat2/)

Install:

```console
$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
$ unzip hisat2-2.1.0-Linux_x86_64.zip
$ ls hisat2-2.1.0/
```

Your directory name may differ, depending on the version of Hisat2 you have.

Update `PATH` and check `hisat2` tool is available:

```console
$ export PATH=~/hisat2-2.1.0:$PATH
$ hisat2 --version
/home/ubuntu/hisat2-2.1.0/hisat2-align-s version 2.1.0
64-bit
Built on login-node03
Wed Jun  7 15:53:42 EDT 2017
Compiler: gcc version 4.8.2 (GCC)
Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
```

Your version of Hisat2 may differ from that shown.

---

## Bowtie

Web site: [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)

**Note:** We are working to add back a Bowtie alignment option.

Install:

```console
$ wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip/download -O bowtie.zip
$ unzip bowtie.zip
$ ls bowtie-1.2.2-linux-x86_64/
```

Your directory name may differ, depending on the version of Bowtie you have.

Update `PATH` and check `bowtie` tool is available:

```console
$ export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
$ bowtie --version
/home/ubuntu/bowtie-1.2.2-linux-x86_64/bowtie-align-s version 1.2.2
64-bit
Built on 462e5beae518
Mon Dec 11 19:27:01 UTC 2017
Compiler: gcc version 4.8.2 20140120 (Red Hat 4.8.2-15) (GCC)
Options: -O3 -m64  -Wl,--hash-style=both -DWITH_TBB
-DPOPCNT_CAPABILITY -g -O2 -fvisibility=hidden -I/hbb_exe/include
-g -O2 -fvisibility=hidden -I/hbb_exe/include
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
```

Your version of Bowtie may differ from that shown.

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

**Note:** Release 2.14.0 or later is required as it includes the [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html) package.

We recommend using R 3.4+.

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

### Check R has installed

```console
$ R --version
R version 3.5.1 (2018-07-02) -- "Feather Spray"
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

### Bioconductor Rsamtools, rdf5, rtracklayer, Biostrings, ShortRead

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

## Check names and versions of Python packages

Run:

```console
$ conda list
```

Alternatively, run:

```console
$ pip list
```

The Python packages and their versions will be listed:

```
biopython                 1.73

cutadapt                  1.18

gffutils                  0.9

h5py                      2.9.0

pysam                     0.15.2

pyyaml                    5.1

samtools                  1.9

umi_tools                 1.0.0
```

Your versions may differ from those shown.

---

## Check names and versions of R packages

(from [list user installed packages](https://www.r-bloggers.com/list-of-user-installed-r-packages-and-their-versions/))

Either run bash script:

```console
$ Rscript install/list-r-packages.R
```

Or run in R:

```r
> ip <- as.data.frame(installed.packages()[,c(1,3:4)])
> rownames(ip) <- NULL
> ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
> print(ip, row.names=FALSE)
```

The R packages and their versions will be listed:

```
            Biostrings    2.46.0

            ggplot2     3.1.1

            optparse     1.6.2

            plotly     4.9.0

            RcppRoll     0.3.0

            readr     1.3.1

            rhdf5    2.22.0

            Rsamtools    1.30.0

            rtracklayer    1.38.3

            shiny     1.3.2

            ShortRead    1.36.1

            tidyr     0.8.3
```

Your versions may differ from those shown.

---

## RiboViz

Get RiboViz:

```console
$ git clone https://github.com/riboviz/riboviz
```

---

## Tested platforms

These instructions were tested on:

| Operating System | Memory (GB) | Processors | RAM (GB) | Python | R |
| ---------------- | ----------- | ---------- | -------- | ------ | - |
| Ubuntu 18.04 | 8 | 4 | 20  | 3.7.3 | 3.4.4 |
| CentOS 7.4.1708 (Core) | 8 | 4 | 20 | 3.7.3 | 3.5.2 |
