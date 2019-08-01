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

Web site: [bedtools](http://bedtools.readthedocs.io/en/latest/)

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

**Check**

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

---

## Python

Web site: [python](https://www.python.org/)

If you already have Python you can skip this step. If you don't have Python then we recommend [Miniconda](https://conda.io/miniconda.html) Python.

Either Python 2.7+ or Python 3.6+ can be used.

### Miniconda Python 2.7

Install:

```console
$ wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda2.sh
$ bash miniconda2.sh -b -p $HOME/miniconda2
```

Activate environment and check:

```console
$ source $HOME/miniconda2/bin/activate
$ python -V
Python 2.7.16 :: Anaconda, Inc.
```

### Miniconda Python 3.6

Install:

```console
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
$ bash miniconda3.sh -b -p $HOME/miniconda3
```

Activate environment and check:

```console
$ source $HOME/miniconda3/bin/activate
$ python -V
Python 3.7.3
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

### Cutadapt

Web sites:

* [GitHub](https://github.com/marcelm/cutadapt)
* [readthedocs](https://cutadapt.readthedocs.io/)

Install:

```console
$ conda install -y -c bioconda cutadapt 
```

Check:

```console
$ conda list | grep cutadapt
cutadapt                  2.3             ...
$ cutadapt --v
2.3
```

**Note:** for Python 2.7 the version could be 1.18. It is OK to use this version.

### pysam

Web sites:

* [GitHub](https://github.com/pysam-developers/pysam/)
* [readthedocs](https://pysam.readthedocs.io/)

Install:

```console
$ conda install -y -c bioconda pysam
```

Check:

```console
$ conda list | grep pysam
pysam                     0.15.2 ...
$ samtools --version
samtools 1.9
Using htslib 1.9
Copyright (C) 2018 Genome Research Ltd.
```

### BioPython

Web site:

* [Biopython](http://biopython.org/)

Install:

```console
$ conda install -y -c anaconda biopython
```

Check:

```console
$ conda list | grep biopython
biopython                 1.73             ...
```

### gffutils

Web site:

* [gffutils](http://daler.github.io/gffutils/)

Install:

```console
$ pip install gffutils
```

Check:

```console
$ conda list | grep gffutils
gffutils                  0.9                      ...
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

Check:

```console
$ conda list | grep h5py
h5py                      2.9.0            ...
$ python
```
```python
>>> import h5py
>>> h5py.run_tests()
...
OK (skipped=16, expected failures=6)
```

---

## Hisat2

Web site: [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)

Install:

```console
$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
$ unzip hisat2-2.1.0-Linux_x86_64.zip 
$ ls hisat2-2.1.0/
```

Update `PATH` and check:

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

Update `PATH` and check:

```console
$ export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
$ bowtie --version
/home/ubuntu/bowtie-1.2.2-linux-x86_64/bowtie-align-s version 1.2.2
64-bit
Built on 462e5beae518
Mon Dec 11 19:27:01 UTC 2017
Compiler: gcc version 4.8.2 20140120 (Red Hat 4.8.2-15) (GCC) 
Options: -O3 -m64  -Wl,--hash-style=both -DWITH_TBB -DPOPCNT_CAPABILITY -g -O2 -fvisibility=hidden -I/hbb_exe/include   -g -O2 -fvisibility=hidden -I/hbb_exe/include  
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
```

## Create `setenv.sh` to configure Hisat2 and Bowtie paths

Create a `setenv.sh` script with contents:

```console
#!/usr/bin/env bash
export PATH=~/hisat2-2.1.0:$PATH
export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
```

---

## R 2.14.0+

Web sites:

* [The R Project for Statistical Computing](https://www.r-project.org/)
* [The Comprehensive R Archive Network](https://cran.r-project.org/) (CRAN).

**Note:** Release 2.14.0 or later is required as it includes the [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html) package.

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

### Check

```console
$ R --version
R version 3.5.1 (2018-07-02) -- "Feather Spray"
```

---

## R packages

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
* Click the link of the specific packag.

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
Error in system2(file.path(R.home("bin"), "R"), c(if (nzchar(arch)) paste0("--arch=",  : 
  cannot popen ' '/usr/lib/R/bin/R' --no-save --slave 2>&1 < '/tmp/Rtmpw3pOH7/file12471113d0d2b'', probable reason 'Cannot allocate memory'
* removing "/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.5/Rsamtools"
Warning in q("no", status = status, runLast = FALSE) :
  system call failed: Cannot allocate memory

The downloaded source packages are in
        "/tmp/RtmpOEVbaL/downloaded_packages"
installation path not writeable, unable to update packages: foreign
Warning message:
In install.packages(pkgs = doing, lib = lib, ...) :
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

The Python packages will be listed:

```
biopython                 1.73

cutadapt                  1.18

gffutils                  0.9

h5py                      2.9.0

pysam                     0.15.2

pyyaml                    5.1
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

The R packages will be listed:

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
$ git clone https://github.com/riboviz/RiboViz
```

---

## Tested platforms

These instructions were tested on:

| Operating System | Resources | Python | R |
| ---------------- | --------- | ------ | - |
| Ubuntu 18.04 | 8GB memory, 4 processors, 20GB RAM | 2.7.16 | 3.4.4 |
| Ubuntu 18.04 | 8GB memory, 4 processors, 20GB RAM | 3.7.3 | 3.4.4 |
| CentOS 7.4.1708 (Core) | 8GB memory, 4 processors, 20GB RAM | 2.7.16 | 3.5.2 |
| CentOS 7.4.1708 (Core) | 8GB memory, 4 processors, 20GB RAM | 3.7.3 | 3.5.2 |
