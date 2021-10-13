# Adding, using, renaming, and removing configuration parameters

* [Adding new configuration parameters](#adding-new-configuration-parameters)
* [Using directory configuration parameters in Nextflow](#using-directory-configuration-parameters-in-Nextflow)
* [Renaming configuration parameters](#renaming-configuration-parameters)
* [Removing configuration parameters](#removing-configuration-parameters)

## Adding new configuration parameters

To add a new configuration parameter, make the following updates.

### 1. Update configuration files

Add the new parameter and its default value to each of the following example configuration files:

```
vignette/remote_vignette_config.yaml
vignette/riboviz_env.yaml
vignette/simdata_multiplex_config.yaml
vignette/simdata_umi_config.yaml
vignette/vignette_config.yaml
riboviz/test/config/vignette_config_current.yaml
```

### 2. Update Python code

Update `riboviz/params.py`, creating a Python constant for the new parameter. For example:

```python
EXAMPLE_CONFIG = "example_config"
""" Example configuration parameter. """
```

Update `riboviz/upgrade_config.py`, adding the new parameter and its default value to the `UPDATES` dictionary. For example:

```python
UPDATES = {
    ...
    params.EXAMPLE_CONFIG = 64,
    ...
}
```

If there is no default value then use the value `None`.

Run the tests for `riboviz/upgrade_config.py` and check that they all pass:

```console
$ pytest riboviz/test/test_upgrade_config.py 
```

### 3. Update Nextflow, `prep_riboviz.nf`

Add the parameter and its description to the `helpMessage()` text.

If the parameter is an optional parameter then, in the block starting with comment `Initialise and validate configuration`, add a line defining the default value for the parameter. For example:

```
params.example_config = 64
```

If the parameter is a mandatory parameter, which **does not** refer to a file or directory path, then add a condition to check that it has been defined and, if not, exit with an error. For example:

```
if (! params.containsKey('example_config')) {
    exit 1, "Undefined adapters (example_config)"
}
```

If the configuration parameter is a mandatory directory parameter, for a directory that must exist before the workflow runs, use the following pattern (replacing `example_dir` with name of the parameter and updating the error message):

```
if ((! params.containsKey('example_dir')) || (! params.example_dir)) {
    exit 1, "Undefined example directory (example_dir)"
}
example_dir = replace_tokens(params.example_dir, riboviz_env_paths)
```

Else if the configuration parameter is an optional directory parameter, use the following pattern (replacing `example_dir` with the name of the parameter):

```
example_dir = replace_tokens(params.example_dir, riboviz_env_paths)
```

Else if the configuration parameter is a mandatory file parameter, use the following pattern (replacing `example_file` with the new parameter name and updating the error message):

```
if ((! params.containsKey('example_file')) || (! params.example_file)) {
    exit 1, "Undefined example TSV file (example_file)"
}
example_file = file(replace_tokens(params.example_file,
                                   riboviz_env_paths))
example_tsv = Channel.fromPath(example_file, checkIfExists: true)
```

Else if the configuration parameter is an optional file parameter, use the following pattern (replacing `example_file` with the new parameter name and updating the error message and the argument to `file`):

```
if (params.containsKey('example_file') && params.example_file) {
    example_file = file(replace_tokens(params.example_file,
                                       riboviz_env_paths))
    if (! example_file.exists()) {
        exit 1, "No such example file (example_file): ${example_file}"
    }
    example_tsv = Channel.fromPath(example_file,
                                   checkIfExists: true)
    is_example_file = true
} else {
    example_tsv = file("missing_example_file")
    is_example_file = false
}
```

### 4. Update documentation

Update `docs/user/prep-riboviz-config.md`, adding an entry for the parameter in the table within the `Configuration parameters` section.

If the new parameter is a file or directory parameter, update `docs/user/prep-riboviz-config.md`, adding the new parameter to the `Configuring file paths and directories` section.

If the new parameter is a directory parameter, update `docs/developer/how-tos.md`, adding the new parameter to the list in the `Using directory configuration parameters in Nextflow` section.

---

## Using directory configuration parameters in Nextflow                        

Some configuration parameters take values that are absolute or relative paths to directories. These are listed in the [Configuring file paths and directories](../user/prep-riboviz-config.md#configuring-file-paths-and-directories) section in [Configuring the riboviz workflow](../user/prep-riboviz-config.md) and the ones that are directories are reproduced below:                                 
* `dir_in`
* `dir_index`
* `dir_tmp`
* `dir_out`

To allow for flexibility as to where users\' input files are located, such configuration variables allow for the, optional, use of environment variable tokens in their values (see [Environment variables and configuration tokens](../user/prep-riboviz-config.md#environment-variables-and-configuration-tokens)).
                                                                              
Token replacement is done as part of the configuration validation code at the top of `prep_riboviz.nf`. Each such parameter, `params.<PARAM>`, is converted into a new variable, `<PARAM>`.

If writing Nextflow code or processes that uses one of these configuration parameters above, then use `<PARAM>` and not `params.<PARAM>` in your new code. For example, if desiring to use the configuration parameter `dir_out` then use `dir_out` and not `params.dir_out`. The processing of `params.dir_out` and any token replacement will already have been done for you by the developer who added the parameter.

---

## Renaming configuration parameters

To rename a configuration parameter, make the following updates.

### 1. Update configuration files

Rename the parameter in each of the following example configuration files:

```
vignette/remote_vignette_config.yaml
vignette/riboviz_env.yaml
vignette/simdata_multiplex_config.yaml
vignette/simdata_umi_config.yaml
vignette/vignette_config.yaml
riboviz/test/config/vignette_config_current.yaml
```

### 2. Update Python code

Update `riboviz/params.py`, renaming the parameter and the associated Python constant for the parameter. For example, if renaming `old_example` to `new_example`, and `riboviz.params` defines:

```python
OLD_EXAMPLE = "old_example"
""" Example configuration parameter. """
```

then `riboviz.params` would be updated to:

```python
NEW_EXAMPLE = "new_example"
""" Example configuration parameter. """
```

Update `riboviz/upgrade_config.py`:

* Add an entry for the old parameter name and its new name to the `Configuration parameters that have been renamed are updated` list in the module comment at the top of the file.
* Add an entry for the old parameter name and its new name in the `RENAMES` dictionary.
* Rename the parameter and to its new name in the `UPDATES` dictionary.

For example, if renaming `old_example` to `new_example`, and `riboviz.upgrade_config` defines:

```python
UPDATES = {
    ...
    params.OLD_EXAMPLE = 64,
    ...
}
```

then `riboviz.upgrade_config` would be updated to:

```python
RENAMES = {
    ...
    "old_example": params.NEW_EXAMPLE,
    ...
}

UPDATES = {
    ...
    params.NEW_EXAMPLE = 64,
    ...
}

```

Run the tests for `riboviz/upgrade_config.py` and check that they all pass:

```console
$ pytest riboviz/test/test_upgrade_config.py 
```

Search for all references to the `riboviz.params` constant as it was prior to renaming, across all the Python code, and update these references.

### 3. Update Nextflow, `prep_riboviz.nf`

Rename the parameter in the `helpMessage()` text.

Rename the parameter and any references to it in the body of `prep_riboviz.nf`.

Rename any Nextflow variables that include the parameter name. For example, if renaming a file configuration parameter from `old_file` to `new_file` and there is a complementary channel called `old_file_tsv` then this would be renamed to `new_file_tsv`.

### 4. Update documentation

Search for all occurrences of the parameter in the documentation and rename these.

### 5. Rerun all tests

* [Run Python tests and workflow tests](./testing.md#run-python-tests-and-workflow-tests).
* [Run R tests](./testing.md#run-r-tests).
* [Run vignette integration tests](./testing.md#run-vignette-integration-tests).

---

## Removing configuration parameters

To remove a configuration parameter that is no longer used, make the following updates.

### 1. Update configuration files

Remove the parameter in each of the following example configuration files:

```
vignette/remote_vignette_config.yaml
vignette/riboviz_env.yaml
vignette/simdata_multiplex_config.yaml
vignette/simdata_umi_config.yaml
vignette/vignette_config.yaml
riboviz/test/config/vignette_config_current.yaml
```

### 2. Update Python code

Update `riboviz/params.py`, removing the parameter and the associated Python constant for the parameter.

Update `riboviz/upgrade_config.py`:

* Remove any entry for the parameter from the `Configuration parameters that have been renamed are updated` list in the module comment at the top of the file.
* Add an entry for the parameter to the `Configuration parameters that are now unused are removed` list in the module comment at the top of the file.
* Remove any occurrence of the parameter from the `UPDATES` and `RENAMES` dictionaries.
* Add an entry for the parameter n the `UNUSED` dictionary.

For example, if removing `deprecated_example`, then `UNUSED` would be updated to:
```python
UNUSED = {
    ...
    "deprecated_example",
    ...
}

```

Run the tests for `riboviz/upgrade_config.py` and check that they all pass:

```console
$ pytest riboviz/test/test_upgrade_config.py 
```

### 3. Update Nextflow, `prep_riboviz.nf`

Remove the parameter from the `helpMessage()` text.

### 4. Update documentation

Search for all occurrences of the parameter in the documentation and remove these.

### 5. Rerun all tests

* [Run Python tests and workflow tests](./testing.md#run-python-tests-and-workflow-tests).
* [Run R tests](./testing.md#run-r-tests).
* [Run vignette integration tests](./testing.md#run-vignette-integration-tests).
