# Deploy a CentOS Virtual Machine using VMWare on Windows 10

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

## Download CentOS 7

Download CentOS 7:

* Visit [CentOS Linux](https://www.centos.org/download/).
* Under CentOS Linux, click 7 (2000).
* Under CentOS Linux, click x84_64.
* Click a mirror URL.
* Click CentOS-7-x86_64-DVD-2009.iso.

The download details are as follows:

* Name: CentOS-7-x86_64-DVD-2009.iso.
* Size: 4712300544 bytes (4494 MB, 4.4 GB)
* SHA256 checksum: e33d7b1ea7a9e2f38c8f693215dd85254c3a4fe446f93f563279715b68d07987  

---

## Create CentOS virtual machine

Create VM:

* Start VMWare Workstation Player.
* Click Create a New Virtual Machine.
* New Virtual Machine Wizard appears.
* Select Installer disc image file (iso).
* Click Browse...
* Use file browser to browse to the `.iso` file e.g. `CentOS-7-x86_64-DVD-2009.iso`.
* Click OK.
* Click Next.
* Name the Virtual Machine appears.
* Enter Virtual machine name: CentOS 7 64-bit
* Select Location: Use the default location. This will be something like `C:\Users\<WINDOWS USER NAME>\<VIRTUAL MACHINE NAME>` e.g. `C:\Users\Jo Smith\Documents\Virtual Machines\CentOS 7 64-bit`.
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

VMWare Workstation Player will now build the CentOS VM:

* `Automatic boot in .. seconds...` appears. Let this count down to 0 and the automatic boot will start.
* `Press the <ENTER> key to begin the installation process.` appears. Wait and the process will begin automatically.
* WELCOME TO CENTOS 7 appears.
* Click preferred language e.g. English (United Kingdom).
* Click Continue.
* INSTALLATION SUMMARY appears.
* Click INSTALLATION DESTINATION.
* INSTALLATION DESTINATION appears.
* Click Done.
* Click SOFTWARE SELECTION.
* SOFTWARE SELECTION appears.
* Select GNOME Desktop on left.
* Select Compatibility Libraries on right.
* Select Development Tools on right.
* Select Security Tools on right.
* Click Done.
* Click NETWORK & HOST NAME.
* NETWORK & HOST NAME appears.
* Click OFF => ON.
* Enter Host name: centos-vm
* Click Apply.
* Click Done.
* Click Begin Installation.

While installation proceeds do the following:

* Click ROOT PASSWORD.
* ROOT PASSWORD appears.
* Enter:
  - Root Password: Enter a password.
  - Confirm: Re-enter the password.
* Click Done.
* Click USER 
* Click User Creation.
* CREATE USER apears.
* Enter:
  - Full name: Enter a name e.g. Jo Smith.
  - User name: Enter a user name e.g. josmith.
  - Password: Enter a password.
  - Confirm password: Re-enter the password.
* Click Done.

Wait for installation to complete. When installation completes:

* Click Reboot.
* Wait for VM to reboot.
* INITIAL SETUP appears.
* Click LICENSE INFORMATION.
* LICENSE INFORMATION appears.
* Check I accept the licence ageement.
* Click Done.
* Click FINISH CONFIGURATION.

Once complete the CentOS prompt will appear. To log in:

* Click user name.
* Enter password.

The desktop will show a Welcome! wizard:

* Select preferred language.
* Click Next.
* Typing appears.
* Select preferred keyboard layout.
* Click Next.
* Privacy appears.
* Set Location Services to OFF.
* Click Next.
* Connect Your Online Accounts appears.
* Click Skip.
* You're ready to go! appears.
* Click Start Using CentOS Linux.

Open a Terminal (bash shell) window:

* Right-click desktop => select Open terminal.

Give yourself `sudo` access so you can have admin priviledges:

```console
$ su
Password: Enter root password.
$ export EDITOR=nano
$ visudo
```

* Scroll down and uncomment:

```
# %wheel  ALL=(ALL)       NOPASSWD: ALL
```

* i.e. change it to:

```
%wheel  ALL=(ALL)       NOPASSWD: ALL
```

* Click CTRL-X to exit.
* Click Y to save updates.
* Run:

```console
$ usermod -aG wheel USER
```

* Click Power button on top-right of desktop.
* Click Power icon.
* Click Restart.
* Log in.
* Open a Terminal (bash shell) window
* Check sudo access:

```console
$ sudo su -
[root@centos-vm ~]# 
```

* Type CTRL-D to exit.

---

## About these instructions

These instructions were tested on a Dell Latitude E7390 laptop with:

* 64-bit Intel Core i5-8350U CPU 1.7GHz, 1.90GHz, 4 cores, 8 logical processors.
* 16 GB RAM.
* 475 GB hard disk.
* Windows 10 Pro.
