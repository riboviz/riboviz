# Upgrade configuration files from RiboViz 1.x

RiboViz has undergone extensive changes from 1.x. This includes the prerequisites RiboViz needs, how RiboViz is used and also the structure of RiboViz's YAML configuration files. Consult the appropriate sections of the documentation for details.

To assist your migration to the current version of RiboViz, we have included a tool to help you upgrade your 1.x-compliant YAML configuration files, `riboviz.tools.upgrade_config_file`.

The tool can be used as follows:

```console
$ python -m riboviz.tools.upgrade_config_file -i <INPUT_FILE> [-o <OUTPUT_FILE>]
```

where:

* `<INPUT_FILE>` is your current YAML configuration file
* `<OUTPUT_FILE>` is a file into which to write the upgraded YAML configuration file. If omitted then the upgraded content is printed to standard out.

For example:

```console
$ python -m riboviz.tools.upgrade_config_file \
    -i /home/user/riboviz-data/config.yaml -o nu_config.yaml 
```

**Tip:** we strongly recommend you visually inspect the updated configuration to ensure the updates reflect your local environment. This especially relates to any file paths.

For details of the changes applied, see [riboviz.tools.upgrade_config_file](../../riboviz/upgrade_config.py).
