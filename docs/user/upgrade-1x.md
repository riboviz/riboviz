# Upgrade configuration files from RiboViz 1.x

RiboViz has undergone extensive changes from 1.x. This includes the prerequisites RiboViz needs, how RiboViz is used and also the structure of RiboViz's YAML configuration files. Consult the appropriate sections of the documentation for details.

To assist your migration to the current version of RiboViz, we have included a tool to help you upgrade your 1.x-compliant YAML configuration files, `upgrade_config_file`.

The tool can be used as follows:

```console
$ upgrade_config_file -i <INPUT_FILE> [-o <OUTPUT_FILE>]
```

where:

* `<INPUT_FILE>` is your current YAML configuration file
* `<OUTPUT_FILE>` is a file into which to write the upgraded YAML configuration file. If omitted then the upgraded content is printed to standard out.

For example:

```console
$ upgrade_config_file -i /home/user/riboviz-data/config.yaml -o nu_config.yaml 
```

**Note:** it is strongly recommended that you visually inspect the updated configuration to ensure the updates reflect your local environment. This especially relates to any file paths.

**Note:** One specific change you may need to make, if your configuration files date from before release 1.1.0, is as follows. If the path of your `features_file` is `scripts/yeast_features.tsv` then this needs to be updated to `data/yeast_features.tsv`, to reflect the relocation of that file in the repository.

For details of the changes applied, see [riboviz.upgrade_config](../../riboviz/upgrade_config.py).
