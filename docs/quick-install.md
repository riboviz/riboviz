# Quick install scripts

`scripts/` contains simple bash scripts to automate some of the installation of **riboviz**'s dependencies. They are the manual commands of [Install prerequisites](./install.md) in bash script form.

**Note:** These scripts are not robust and do not perform any error handling.

These instructions were written for and tested on Ubuntu 18.04 LTS bionic.

Installing some of these tools requires you to have `sudo` access to install and configure software (or a local system administrator to do this for you).

Authenticate with `sudo`:

```
sudo su -
CTRL-D
```

Run the installer scripts:

```
source install.sh
source install-py.sh
source install-riboviz.sh
source install-r.sh
```

Start R:

```
R
```

Run the R installer script:

```
source("install-r.R")
```

When prompted:

```
Would you like to use a personal library instead? (yes/No/cancel) 
```

Enter:

```
yes
```

When complete, check that R's library paths include your personal library:

```
.libPaths()
```

You should see something like the following (with `/home/ubuntu/` replaced by your home directory):

```
/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.5
```
