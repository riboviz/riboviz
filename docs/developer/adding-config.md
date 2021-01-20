# Adding configuration parameters

If a new configuration parameter is added to the YAML configuration file then the following updates need to be done.

1. Add an entry for the parameter in the "Configuration parameters" table in [Configuring the RiboViz workflow](.../user/prep-riboviz-config.md).
2. Add an entry for the parameter, and its default value, to each of:

```
vignette/vignette_config.yaml
vignette/simdata_umi_config.yaml
vignette/simdata_multiplex_config.yaml
riboviz/test/config/vignette_config_current.yaml
riboviz/test/config/simdata_umi_config_current.yaml
riboviz/test/config/simdata_multiplex_config_current.yaml
```

3. Add a Python constant for the parameter to `riboviz/params.py`:

For example, if the parameter is `some_setting` then add:

```python
SOME_SETTING = "some_setting"
```

4. Update the "Upgrade configuration files" code, `riboviz/upgrade_config.py`:

Add the parameter and its default value to the list under the comment:

```
Expected parameters added between release 2.0 and current are added
along with default values, if they are not already present in the
configuration:
```

Add the parameter to the `UPDATES_20_CURRENT` dictionary. For example, if the parameter was `some_setting` as above then add a mapping from its associated `riboviz.params` value to its default value. For example, if the default value is "some_value", add:

```
params.SOME_SETTING: "some_value"
```

5. Run the tests for `riboviz/upgrade_config.py` and ensure they all pass:

```console
$ pytest riboviz/test/test_upgrade_config.py 
```

These use the YAML configuration files in `vignette/` and `riboviz/test/config/`.
