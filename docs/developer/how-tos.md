# Developer how-tos

* [Clone conda environments](#clone-conda-environments)
* [Understanding YAML `NULL` and Python `None`](#understanding-yaml-null-and-python-none)

---

## Clone conda environments

If you want to install a new Python package to explore it without changing your current conda environment you can first clone your conda environment, then install the package into the clone. For example:

```console
$ conda create --name riboviz-test-install --clone riboviz
$ conda activate riboviz-test-install
```

To remove an environment when you are finished with it you can run, for example:

```console
$ conda deactivate riboviz-test-install
$ conda env remove --name riboviz-test-install
```

For more information, see conda's [Managing environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

---             

## Understanding YAML `NULL` and Python `None`

If a parameter in a YAML file has value `null`, `NULL` or no value at all then, after reading the file into Python (using the `yaml` library), it will have value `None`. For example, given a YAML file with:

```yaml
a: null
b: NULL
c:
```

Loading this with

```python
with open("config.yml") as f:
    config = yaml.load(f, Loader=yaml.SafeLoader)
```

would result in the following all having value `None`:

```python
config['a']
config['b']
config['c']
```
