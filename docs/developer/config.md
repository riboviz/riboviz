# Using and adding configuration parameters

## Using directory configuration parameters in Nextflow

Some configuration parameters take values that are absolute or relative paths to directories. These are listed in the [Configuring file paths and directories](../user/prep-riboviz-config.md#configuring-file-paths-and-directories) section in [Configuring the RiboViz workflow](../user/prep-riboviz-config.md) and the ones that are directories are reproduced below:

* `dir_in`
* `dir_index`
* `dir_tmp`
* `dir_out`

To allow for flexibility as to where users' input files are located, such configuration variables allow for the, optional, use of environment variable tokens in their values (see [Environment variables and configuration tokens](../user/prep-riboviz-config.md#environment-variables-and-configuration-tokens)).

Token replacement is done as part of the configuration validation code at the top of `prep_riboviz.nf`. Each such parameter, `params.<PARAM>`, is copied to a new variable `<PARAM>`.

If writing Nextflow code or processes that uses one of these configuration parameters above, then use `<PARAM>` and not `params.<PARAM>` in your new code. For example, if desiring to use the configuration parameter `dir_out` then use `dir_out` and not `params.dir_out`. The processing of `params.dir_out` and any token replacement will already have been done for you.

---

## Adding configuration parameters

If a new configuration parameter is added to the YAML configuration file and to the Nextflow workflow then the following updates need to be done.

1. Add a Python constant for the parameter to `riboviz/params.py`. For example:

```python
EXAMPLE_FILE = "example_file"
""" Example file. """
EXAMPLE_THRESHOLD = "example_threshold"
""" Example threshold. """
```

2. Add an entry for the parameter, and its default value, to each of:

```
vignette/vignette_config.yaml
vignette/simdata_umi_config.yaml
vignette/simdata_multiplex_config.yaml
vignette/remote_vignette_config.yaml
riboviz/test/config/vignette_config_current.yaml
riboviz/test/config/simdata_umi_config_current.yaml
riboviz/test/config/simdata_multiplex_config_current.yaml
```

3. Update the "Upgrade configuration files" code, `riboviz/upgrade_config.py`:

Add the parameter and its default value to the module comment. If there is no default value then use the value `null`. For example:

```
* ``example_file: null``
* ``example_threshold: 64``
```

Add the parameter to the `UPDATES` dictionary. For example, if the parameter was `some_setting` as above then add a mapping from its associated `riboviz.params` value to its default value. If there is no default value then use the value `None`. For example:

```python
params.EXAMPLE_FILE: None,
params.EXAMPLE_THRESHOLD: 64,
```

4. Run the tests for `riboviz/upgrade_config.py` and ensure they all pass:

```console
$ pytest riboviz/test/test_upgrade_config.py 
```

These use the YAML configuration files in `vignette/` and `riboviz/test/config/`.

5. Add an entry for the parameter in the [Configuration parameters](../user/prep-riboviz-config.md#configuration-parameters) table in [Configuring the RiboViz workflow](../user/prep-riboviz-config.md).

6. If the parameter is optional and has a default value then add the parameter and its default value to `prep_riboviz.nf` e.g.

```
params.example_threshold = 64
```

7. If the configuration parameter takes a value that is an absolute or relative path to a file or directory then see [Adding new path configuration parameters in Nextflow](#adding-new-path-configuration-parameters-in-nextflow) below.

---

## Adding path configuration parameters in Nextflow

As described in [Using directory configuration parameters in Nextflow](#using-directory-configuration-parameters-in-nextflow), a number of configuration parameters take values that are absolute or relative paths to files or directories. Token replacement is done as part of the configuration validation code at the top of `prep_riboviz.nf`. Adding Nextflow code to handle new path parameters can be done as follows.

For a mandatory directory parameter with a directory that must exist before the workflow runs, use the following pattern (replacing `dir_in` with the new parameter name and updating the error message):

```
if (! params.containsKey('dir_in')) {
    exit 1, "Undefined input directory (dir_in)"
}
dir_in = replace_tokens(params.dir_in, riboviz_env_paths)
```

For an optional directory parameter, use the following pattern (replacing `dir_out` with the new parameter name):

```
dir_out = replace_tokens(params.dir_out, riboviz_env_paths)
```

For a mandatory file parameter, use the following pattern (replacing `rrna_fasta_file` with the new parameter name and updating the error message):

```
if (! params.containsKey('rrna_fasta_file')) {
    exit 1, "Undefined rRNA FASTA file (rrna_fasta_file)"
}
rrna_fasta_file = file(replace_tokens(params.rrna_fasta_file,
                                      riboviz_env_paths))
rrna_fasta = Channel.fromPath(rrna_fasta_file, checkIfExists: true)
```

For an optional file parameter, use the following pattern (replacing `asite_disp_length_file` with the new parameter name and updating the error message):

```
if (params.containsKey('asite_disp_length_file')
    && params.asite_disp_length_file) {
    asite_disp_length_file = file(replace_tokens(params.asite_disp_length_file,
                                                 riboviz_env_paths))
    if (! asite_disp_length_file.exists()) {
        exit 1, "No such A-site displacement file (asite_disp_length_file): ${asite_disp_length_file}"
    }
    asite_disp_length_txt = Channel.fromPath(asite_disp_length_file,
                                             checkIfExists: true)
    is_asite_disp_length_file = true
} else {
    asite_disp_length_txt = file("Missing_aside_disp_length_file")
    is_asite_disp_length_file = false
}
```

The new parameter should also be added to the list in the [Configuring file paths and directories](../user/prep-riboviz-config.md#configuring-file-paths-and-directories) section in [Configuring the RiboViz workflow](../user/prep-riboviz-config.md).

If the parameter is a directory then add it to the list in [Using directory configuration parameters in Nextflow](#using-directory-configuration-parameters-in-nextflow) above.
