# Developer how-tos

## Clone conda environments

If you want to install a new Python package to explore it without changing your current conda environment you can first clone your cona environment, then install the package into the clone. For example:

```console
$ conda create --name riboviz-test-install --clone riboviz
$ conda activate riboviz-test-install
```

---

## Use directory configuration parameters in Nextflow                        

Some configuration parameters take values that are absolute or relative paths to directories. These are listed in the [Configuring file paths and directories](../user/prep-riboviz-config.md#configuring-file-paths-and-directories) section in [Configuring the riboviz workflow](../user/prep-riboviz-config.md) and the ones that are directories are reproduced below:                                 
* `dir_in`
* `dir_index`
* `dir_tmp`
* `dir_out`

To allow for flexibility as to where users' input files are located, such configuration variables allow for the, optional, use of environment variable tokens in their values (see [Environment variables and configuration tokens](../user/prep-riboviz-config.md#environment-variables-and-configuration-tokens)).
                                                                              
Token replacement is done as part of the configuration validation code at the top of `prep_riboviz.nf`. Each such parameter, `params.<PARAM>`, is converted into a new variable, `<PARAM>`.

If writing Nextflow code or processes that uses one of these configuration parameters above, then use `<PARAM>` and not `params.<PARAM>` in your new code. For example, if desiring to use the configuration parameter `dir_out` then use `dir_out` and not `params.dir_out`. The processing of `params.dir_out` and any token replacement will already have been done for you by the developer who added the parameter.

---             

## YAML `NULL` and Python `None`

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
