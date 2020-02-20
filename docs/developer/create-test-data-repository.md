# Creating a regression test data repository

## Caution

Any repository hosted on GitHub should not exceed 1GB in size. GitHub's [What is my disk quota?](https://help.github.com/en/github/managing-large-files/what-is-my-disk-quota) comments "If your repository exceeds 1GB, you might receive a polite email from GitHub Support requesting that you reduce the size of the repository to bring it back down."

## Create a repository on GitHub

Create a repository, `test-data-<YYYYMMDD>` within the [riboviz](https://github.com/riboviz) project on GitHub. For example: ``test-data-20200220`.

## Clone a local copy of the repository

```console
$ git clone https://github.com/riboviz/<REPOSITORY>
```

## Copy vignette data into the repository

```console
$ cd riboviz
$ cp vignette/vignette_config.yaml ~/<REPOSITORY>
$ cp -r vignette/logs ~/<REPOSITORY>
$ cp -r vignette/output ~/<REPOSITORY>
```

## Document execution environment

```console
$ source install/environment.sh > ~/<REPOSITORY>/environment.txt 2>&1
```

## Create `README.md`

Copy template `README.md`:

```console
$ cp docs/developer/test-readme-template.md ~/<REPOSITORY>/README.md
```

Edit `<REPOSITORY>/README.md` and fill in the details about the regression test data.

## Add, commit and push data

```console
$ cd ~/<REPOSITORY>
$ git add .
$ git commit -m "Added test data"
$ git push origin master
```
