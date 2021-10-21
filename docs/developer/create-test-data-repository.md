# Creating a test data repository

The process for creating a test data repository is as follows:

1. [Create a repository on GitHub](#1-create-a-repository-on-github)
2. [Clone a local copy of the new repository](#2-clone-a-local-copy-of-the-new-repository)
3. [Copy vignette data into the repository](#3-copy-vignette-data-into-the-repository)
4. [Create `README.md`](#4-create-readme.md)
5. [Add, commit and push data](#5-add-commit-and-push-data)


---

## Caution

Any repository hosted on GitHub should not exceed 1GB in size. GitHub's [What is my disk quota?](https://help.github.com/en/github/managing-large-files/what-is-my-disk-quota) comments "If your repository exceeds 1GB, you might receive a polite email from GitHub Support requesting that you reduce the size of the repository to bring it back down."

---

## 1. Create a repository on GitHub

Create a repository:

* Visit https://github.com/riboviz
* Click New
* Enter Repository name:
  - For data corresponding to a release, enter `test-data-<TAG>` e.g. `test-data-2.1`.
  - For data produced on a given date, enter `test-data-<YYYYMMDD>` e.g. `test-data-20210618`.
* Select Public.
* Select Add a README file.
* Click Create repository.

---

## 2. Clone a local copy of the new repository

```console
$ git clone https://github.com/riboviz/<REPOSITORY>
```

---

## 3. Copy vignette data into the repository

```console
$ cd riboviz
$ cp vignette/vignette_config.yaml ~/<REPOSITORY>
$ cp -r vignette/output ~/<REPOSITORY>
```

---

## 4. Create `README.md`

Copy template `README.md`:

```console
$ cp docs/developer/test-data-readme.md ~/<REPOSITORY>/README.md
```

Edit `<REPOSITORY>/README.md` and fill in the details about the test data.

Document environment under which test data was produced:

```console
$ source bash/environment-tables.sh >> ~/<REPOSITORY>/README.md
```

---

## 5. Add, commit and push data

```console
$ cd ~/<REPOSITORY>
$ git add .
$ git commit -m "Added test data"
$ git push origin main
```
