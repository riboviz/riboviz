# Using the Linux `screen` command

Linux's `screen` command provides a virtual terminal multiplexer. It allows us to run a number of different sessions (or windows, or virtual terminals) withina single terminal, or console, window. This can be useful if we do not have access to a graphical user interface with multiple windows, which may be the case when using high-performance computing systems, for example.

One use case is if we want to have one or more files open in editors while also running programs that use those files - for example editing source code and running a compiler, or editing configuration files and running an analysis program - without having to repeatedly open and close the files within the editors.

This document gives a short introduction to `screen`.

## Install `screen`

To see if `screen` is available on your system, run:

```console
$ screen -v
Screen version 4.01.00devel (GNU) 2-May-06
```

If it is not, if you have permission to run `sudo`, then you can install `screen` as follows. If you don't have `sudo` access you will have to ask a local system administrator to run these commands for you.

* Ubuntu users:

```console
$ apt-get -y install screen
```

* CentOS users:

```console
$ yum -y install screen
```

## A counter program

Create a file, `counter.sh`, that is bash script that prints values up to a given bound:

```
#!/bin/bash
BOUND=$1
for (( i = 0; i < $BOUND; i++ ))
do
    echo "$i"
done
```

Run `counter.sh`:

```console
$ bash counter.sh 1000000
0
1
2
3
...
```

Press CTRL+C to stop `counter.sh`.

## Start, detach, reattach, stop screen sessions

Start a new screen session:

```console
$ screen
```

The commands we now run will be within the context of the screen session we have created. Start `counter.sh`:

```console
$ bash counter.sh 1000000
0
1
2
3
...
```

To detach from current screen session leaving its content running, in this case `counter.sh`, press CTRL+A then press d.

A `detached` message will be shown with the identifier of the screen session that we detached:

```console
[detached from 17476.pts-4.centos]
```

To see a list of all active screens sessions, run:

```console
$ screen -ls
There is a screen on:
	17476.pts-4.centos	(Detached)
1 Socket in /var/run/screen/S-centos.
```

To reattach to a specific screen session, run the following, providing the identifier of the screen:

```console
$ screen -r 17476
...
999997
999998
999999
$
```

If we had left a command running in the screen session, before we detached, as we did in this example, then we will see its current state. For `counter.sh` it will likely have completed, as shown above.

To kill a screen session we are in, we can do one of:

* Press CTRL+A then k. A `Really kill this window [y/n]` prompt will appear.
* Press CTRL+D.
* Run:

```console
$ exit
```

A `[screen is terminating]` message will be shown.

To give a screen session a more memorable name when, we can use `-S`. Run:

```console
$ screen -S counter
$ bash counter.sh 1000000
```

Press CTRL+A then press d to detach. Now list the screens:

```console
$ screen -ls
There is a screen on:
        18107.counter       (Attached)
1 Socket in /var/run/screen/S-centos.
```

To reattach to the screen, run the following, providing the name of the screen:

```console
$ screen -r counter
```

## Use screen session windows

`screen` can support multiple processes within a session. This allows us to run multiple commands concurrently and view their progress. For example we could have a program file open in one window and run our program in another.

**Note:** In this context, a screen session "window" is not a graphical user interface-style window, but rather a view within a screen session. Screen session windows appear within the bounds of a single console or terminal window.

Open `counter.sh` in the `nano` editor:

```console
$ nano counter.sh
```

To create a new screen session window, press CTRL+A then press c. The screen with `nano` should disappear and a new window with a prompt will appear:

```console
$
```

Run `counter.sh`:

```console
$ bash counter.sh 5
0
1
2
3
4
```

To switch back to the window with `nano`, press CTRL+A then press 0 (windows are indexed from 0).

Using `nano`, change the line:

```
    echo "$i"
```

to:

```
    echo "Count: $i"
```

Press CTRL+O to save the change to `counter.sh` (CTRL+O is a `nano` shortcut).

To switch back to the window in which we ran `counter.sh`, press CTRL+A then press 1.

Rerun `counter.sh`:

```console
$ bash counter.sh 5
Count: 0
Count: 1
Count: 2
Count: 3
Count: 4
```

Using CTRL+A and a window name allows us to switch between both windows.

Press CTRL+A then press " to list all windows:

```
 Num Name                                                         Flags

   0 centos@centos:~                                                  $
   1 centos@centos:~                                                  $
```

Each `Name` is a default window name. Yours may differ from those above.

Use the up and down arrows to select a window and then press ENTER to select that window.

Select window 0, the one with `nano`.

### Rename screen session windows

Screen session windows can be given more meaningful names.

Press CTRL+A then press A to rename the current window. You will be prompted for a new name, with the default name shown:

```
Set window's title to: centos@centos:~   
```

Delete the default name and enter a more meaningful one:

```
Set window's title to: Edit counter
```

Press CTRL+A then press 1 to switch to the other window.

Press CTRL+A then press A and rename the current window to `Run counter`

```
Set window's title to: Run counter
```

Press CTRL+A then press " to list all the windows:

```
 Num Name                                                         Flags

   0 Edit counter                                                     $
   1 Run counter                                                      $
```

Now the names are more usable.

Select window 0, `Edit counter`, the one with `nano`.

### Detatch and reattach screen session

Press CTRL+A then press d to detach from the current screen session.

Press CTRL+A and a window number. It will now have no effect as we are no longer in the context of our screen session. But do not worry, our windows are still live within the context of the screen session.

Reattach the screen session:

```console
$ screen -r counter
```

Press CTRL+A and a window number. Our windows are available once more.

### View multiple screen session windows concurrently

We can view multiple screen windows at the same time.

Press CTRL+A, then SHIFT+S (hold down SHIFT, press S). The terminal window will split into two regions. The top region will display the content of the current screen window session, labelled with its number and name `0 Edit counter`. The bottom region will be empty, labelled with `--`.

To switch to the bottom window, press CTRL+A then press TAB.

Press CTRL+A then press 1 to display the window in which we ran `counter.sh`. Alternatively, press CTRL+A then press " and select window 1, `Run Counter`.

Rerun `counter.sh`.

Press CTRL+A then press TAB to switch to the top region.

Change the `echo` message from `"Count: $i"` to `"Progress: $i"`.

Press CTRL+O to save `counter.sh`

Press CTRL+A then press TAB to switch to the bottom region.

Rerun `counter.sh`:

```console
$ bash counter.sh 5
Progress: 0
Progress: 1
Progress: 2
Progress: 3
Progress: 4
```

To exit from the split screen display, press CTRL+A then press Q to close all regions but the current one.

---

## Use screen with SSH

If you SSH into a server you can create a screen session in which to do your work. If you detach from that screen session then exit from SSH, you can log in at a later date and reattach your screen session. This can be useful if running long-running tasks. For example:

```console
$ ssh someuser@eddie.ecdf.ed.ac.uk
$ screen -S sshexample
$ top
```

`top` shows processes running on the current system.

Press CTRL+A then press d to detach the screen session.

```
[detached from 109726.sshexample]
```

Exit from the server:

```console
$ exit
logout
Connection to eddie.ecdf.ed.ac.uk closed.
```

Log in again:

```console
$ ssh someuser@eddie.ecdf.ed.ac.uk
$ screen -ls
There is a screen on:
        109726.sshexample       (Detached)
1 Socket in /var/run/screen/S-mjj.
```

Reattach the detached screen session:

```console
$ screen -r sshexample
```

You should see `top` again.

Press q to exit `top`.

Exit screen session then log out of the server:

```console
$ exit
[screen is terminating]
$ exit
logout
Connection to eddie.ecdf.ed.ac.uk closed.
```

## Addtional screen window commands

To create two regions split vertically, press CTRL+A then press | (pipe).

Regions can be nested. Within a region you can press CTRL+A and S or CTRL+A and | to split that region horizontally or vertically.

To move back and forth between the current and previous window, press CTRL+A then press A.

## CTRL+A gotcha

You will notice that the `screen` commands used start with CTRL+A. This means that other commands that use CTRL+A will not work within the context of a screen session. This includes using CTRL+A to go to the start of a line within a bash shell or within an Emacs editor session.

## Acknowledgements and further information

This document was derived from:

* [How to Use Linux's screen Command](https://www.howtogeek.com/662422/how-to-use-linuxs-screen-command/), Dave McKay, How-To Geek, 27 March, 2020.
* [How To Use Linux Screen](https://linuxize.com/post/how-to-use-linux-screen/), Linuxize, 13 April 2021.

See also GNU [Screen User's Manual](https://www.gnu.org/software/screen/manual/screen.html).
