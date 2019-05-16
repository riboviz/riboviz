# Quick install scripts

`scripts/` contains simple bash scripts to automate some of the installation of **riboviz**'s dependencies. They are the manual commands of [Install prerequisites](./install.md) in bash script form.

These scripts were written for Ubuntu 18.04 and CentOS 7.4.

Running these scripts requires you to have permission to run `sudo` to install and configure software. If you don't have `sudo` access you will have to ask a local system administrator to run these commands for you.

**Note:** These scripts are not robust and do not perform any error handling.

Authenticate with `sudo`:

```
sudo su -
CTRL-D
```

Install operating system packages:

* Ubuntu:

```bash
source install-ubuntu.sh
```

* CentOS:

```bash
source install-centos.sh
```

Install Python:

* **Note**: This script installs both Python 2 and 3. If you only want to install one or the other then edit the file and comment out the code that installs the version you do not want.

```bash
source install-py.sh
```

Install Hisat2 and Bowtie:

```bash
source install-hisat-bowtie.sh
```

Install R:

* Ubuntu:

```bash
source install-r-ubuntu.sh
```

* CentOS

```bash
source install-r-centos.sh
```

Install R packages:

```bash
Rscript install-R.r
```

Install R Bioconductor packages:

* R 3.4 users:

```bash
Rscript install-r-3.4-bioconductor.R
```

* R 3.5 users:

```bash
Rscript install-r-3.5-bioconductor.R
```

When complete, check that R's library paths include your personal library:

```bash
Rscript -e ".libPaths()"
```

You should see something like:

* Ubuntu:

```
[1] "/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.4"
[2] "/usr/local/lib/R/site-library"                 
[3] "/usr/lib/R/site-library"                       
[4] "/usr/lib/R/library"     
```

* CentOS:

```
[1] "/home/centos/R/x86_64-redhat-linux-gnu-library/3.5"
[2] "/usr/lib64/R/library"                              
[3] "/usr/share/R/library"  
```
