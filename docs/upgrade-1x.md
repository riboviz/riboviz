# Upgrade configuration files from RiboViz 1.x

RiboViz has undergone extensive changes from 1.x. This includes the prerequisites RiboViz needs, how RiboViz is used and also the structure of RiboViz's YAML configuration files. Consult the appropriate sections of the documentation for details.

To assist your migration to the current version of RiboViz, we have included a tool to help you upgrade your 1.x-compliant YAML configuration files, `upgrade_config_file.py`.

The tool can be used as follows:

```console
$ python -m riboviz.tools.upgrade_config_file -i <INPUT> [-o <OUTPUT>]
```

where:

* `<INPUT>` is your current YAML configuration file
* `<OUTPUT>` is a file into which to write the upgraded YAML configuration file. If omitted then the upgraded content is printed to standard out.

For example:

```console
$ python -m riboviz.tools.upgrade_config_file \
    -i /home/user/riboviz-data/config.yaml
    -o nu_config.yaml 
```

**Tip:** we strongly recommend you visually inspect the updated configuration to ensure the updates reflect your local environment. This especially relates to any file paths.

## Changes applied by `upgrade_config_file.py`

Configuration parameters that have been renamed from 1.x are updated:

* `rRNA_fasta` => `rrna_fasta_file`
* `orf_fasta` => `orf_fasta_file`
* `rRNA_index` => `rrna_index_prefix`
* `orf_index` => `orf_index_prefix`
* `nprocesses` => `num_processes`
* `MinReadLen` => `min_read_length`
* `MaxReadLen` => `max_read_length`
* `Buffer` => `buffer`
* `PrimaryID` => `primary_id`
* `SecondID` => `secondary_id`
* `StopInCDS` => `stop_in_cds`
* `isTestRun` => `is_test_run`
* `ribovizGFF` => `is_riboviz_gff`
* `t_rna` => `t_rna_file`
* `codon_pos` => `codon_positions_file`

Expected parameters added between release 1.0.0 and 1.1.0 are added along with default values, if they are not already present in the configuration file:

* `do_pos_sp_nt_freq: true`
* `features_file: data/yeast_features.tsv`

Expected parameters added between release 1.1.0 and the current release are added along with default values, if they are not already present in the configuration file:

* `dir_logs: vignette/logs`
* `cmd_file: run_riboviz_vignette.sh`
* `workflow_files_log_file: workflow_files_vignette.tsv`
* `t_rna_file: data/yeast_tRNAs.tsv`
* `codon_positions_file: data/yeast_codon_pos_i200.RData`
* `count_threshold: 64`
* `asite_disp_length_file: data/yeast_standard_asite_disp_length.txt`

The values of parameters `rrna_index_prefix` and `orf_index_prefix` are updated to be file names only, as, these are now assumed to be relative to `<dir_index>`. For example the configuration parameters:

```yaml
rRNA_index: vignette/index/yeast_rRNA
orf_index: vignette/index/YAL_CDS_w_250
```

are updated to:

```yaml
rRNA_index_prefix: yeast_rRNA
orf_index_prefix: YAL_CDS_w_250
```

The value of parameter `features_file` is changed to reflect the relocation of this file in a `scripts/` directory to its new location in a `data/` directory. For example, the configuration parameter:

```yaml
features_file: scripts/yeast_features.tsv
```

is updated to:

```yaml
features_file: data/yeast_features.tsv
```

As another example, the configuration parameter:

```yaml
features_file: /home/user/riboviz/scripts/yeast_features.tsv
```

is updated to:

```yaml
features_file: /home/user/riboviz/data/yeast_features.tsv
```
