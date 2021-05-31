# Deploy a Ubuntu Virtual Machine using VMWare on Windows 10

## Requirements

You will need at least 30GB of free disk space.

---

## Install VMWare Workstation Player

Download VMware Workstation 15.5.6 Player for Windows 64-bit Operating Systems from Download Product,[15.5.6](https://my.vmware.com/en/web/vmware/downloads/details?downloadGroup=PLAYER-1556&productId=800).

The download details are as follows:

* Name: VMware-player-15.5.6-16341506.exe
* Size 152747144 bytes (146 MB)
* SHA256 checksum: 96a8b2da596ec5057f53992b200b834ecba1f6d0ee2290bec2b3e28459c42f7e
* For full details, click Read More on the Download Product page above.

Double-click `VMware-player-15.5.6-16341506.exe` and install VMware Workstation Player.

Notes:

* A newer version, [16.1.2](https://my.vmware.com/web/vmware/downloads/details?downloadGroup=WKST-PLAYER-1612&productId=1039&rPId=66621) (18/05/21), is also available but this has not been tried.

---

## Download Ubuntu

Download Ubuntu 18.04.5 LTS (Bionic Beaver) 64-bit PC (AMD64) desktop image from [Ubuntu 18.04.5 LTS (Bionic Beaver)](https://releases.ubuntu.com/bionic/).

The download details are as follows:

* Name: ubuntu-18.04.5-desktop-amd64.iso
* Size 2193522688 bytes (2092 MB, 2 GB)
* SHA256 checksum: f295570badb09a606d97ddfc3421d7bf210b4a81c07ba81e9c040eda6ddea6a0

Notes:

* The latest versions of Ubuntu are Ubuntu 20.04.2.0 LTS and Ubuntu 21.04, available from [Download Ubuntu Desktop](https://www.ubuntu.com/download/desktop) but these have not been tried.
* Other, older, [releases](https://releases.ubuntu.com/) are also available.

-----

## Create Ubuntu virtual machine

Create VM:

* Start VMWare Workstation Player.
* Click Create a New Virtual Machine.
* New Virtual Machine Wizard appears.
* Select Installer disc image file (iso).
* Click Browse...
* Use file browser to browse to the `.iso` file e.g. `ubuntu-18.04.5-desktop-amd64.iso`.
* Click OK.
* Click Next.
* Easy Install Information appears.
* Enter:
  - Full name: Enter a name e.g. Jo Smith.
  - User name: Enter a user name e.g. josmith.
  - Password: Enter a password.
  - Confirm: Re-enter the password.
* Click Next.
* Name the Virtual Machine appears.
* Enter Virtual machine name: Ubuntu 64-bit 18.04.5
* Select Location: Use the default location. This will be something like `C:\Users\<WINDOWS USER NAME>\<VIRTUAL MACHINE NAME>` e.g. `C:\Users\Jo Smith\Documents\Virtual Machines\Ubuntu 64-bit 18.04.5`.
* Click Next.
* Specify Disk Capacity appears.
* Enter Maximum disk size: 25
* Select Split virtual disk into multiple files.
* Click Next.
* Ready to Create Virtual Machine appears.
* Click Customize Hardware....
* Select Memory for this virtual machine: 8 GB
* Click Processors.
* Number of processor cores: Enter number of processors less than that of your host Windows 10 machine.
* Click Close.
* Click Finish.

VMWare Workstation Player will now build the Ubuntu VM. This many take many minutes. Once complete the Ubuntu prompt will appear. To log in:

* Click user name.
* Enter password.

Do not resize VMWare Workstation Player!

The desktop will show a What's new in Ubuntu wizard:

* Click Next.
* Livepatch appears.
* Click Next.
* Help improve Ubuntu appears.
* Select No, don't send system info.
* Click Next.
* You're ready to go! appears.
* Click Done.

---

## Configure Ubuntu virtual machine

Add Terminal (bash shell) icon to icon bar:

* Right-click desktop => select Open terminal.
* In icon bar on left of desktop, right-click Terminal icon => select Add top Favorites.

Within a terminal window, install Git (enter your password when prompted):

```console
$ sudo apt install -y git
```

Set UK keyboard (optional):

* Click top-right 'V' menu icon, on right-hand side of power icon.
* Click settings (crossed-spanners) button.
* Click Region & Language.
* Click '+'.
* Click English (United Kingdom).
* Click English (UK).
* Click Add.
* Select English (US).
* Click '-'.
* Click 'x' at top right of Settings window to exit.

---

## About these instructions

These instructions were tested on a Dell Latitude E7390 laptop with:

* 64-bit Intel Core i5-8350U CPU 1.7GHz, 1.90GHz, 4 cores, 8 logical processors.
* 16 GB RAM.
* 475 GB hard disk.
* Windows 10 Pro.
