# Developing Nextflow workflow

* [Nextflow style](#nextflow-style)
  - [Process names](#process-names)
  - [Channel names](#channel-names)
  - [Variable names](#variable-names)
* [Adding, updating or removing steps](#adding-updating-or-removing-steps)
* [Debugging](#debugging)
  - [Debugging within Nextflow's `work/` step-specific directories](#debugging-within-nextflows-work-step-specific-directories)
  - [Debugging within a directory emulating Nextflow's `work/` step-specific directories](#debugging-within-a-directory-emulating-nextflows-work-step-specific-directories)

---

## Nextflow style

See also the riboviz [Style guide](./style-guide.md).

### Process names

Nextflow process names must be in CamelCase with the first letter being lower-case. Upper-case is peromitted for acronyms e.g., `ORF`, `CDS`, `APE`, `TPMs`, `RNA`.

For example, `buildIndicesrRNA`, `demultiplex`, `staticHTML`.

### Channel names

Nextflow channel names must be in snake-case i.e., lower-case and delimited by underscores, not hyphens.

For channels which are file names, ensure that the channel name includes the file type as its last component, delimited by an underscore.

For example, `multiplex_sample_sheet_tsv`, `sample_fq`

### Variable names

Nextflow variable names must be in snake-case i.e., lower-case and delimited by underscores, not hyphens.

---

## Adding, updating or removing steps

If adding steps (processes) to the workflow, editing existing steps or removing steps from the workflow then:

* Update the workflow descriptions in `docs/user/prep-riboviz-config.md` to reflect the changes you have made.
* If you have changed configuration parameters, then see [Adding, using, renaming, and removing configuration parameters](./config.md)
* If you have changed temporary or output files, then see [Adding, renaming, and removing temporary or output files](./output-files.md)
* Update the workflow SVG images, `docs/images/*.svg` and see [Updating workflow images](./documentation.md#updating-workflow-images) in [Writing and updating documentation](./documentation.md).

---

## Debugging

Information on Debugging and Nextflow is provided in [Debugging](../user/prep-riboviz-run-nextflow.md#debugging) in [Running the riboviz Nextflow workflow](../user/prep-riboviz-run-nextflow.md). This section contains additional developer-specific information.

### Debugging within Nextflow's `work/` step-specific directories

`.command.sh` scripts within step-specific directories can be rerun within those directories. For example:

```console
$ nextflow run -ansi-log false prep_riboviz.nf -params-file vignette/vignette_config.yaml
N E X T F L O W  ~  version 20.04.1
Launching `prep_riboviz.nf` [big_majorana] - revision: 6c6670470d
...
[38/784d89] Submitted process > collateTpms (WTnone, WT3AT)
Error executing process > 'collateTpms (WTnone, WT3AT)'
...
Work dir:
  /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73
...
$ cd /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73/
$ cat .command.sh 
#!/bin/bash -ue
Rscript --vanilla /home/ubuntu/riboviz/rscripts/collate_tpms.R             --tpms-file=TPMs_all_CDS_all_samples.tsv             WT3AT WT3AT_tpms.tsv WTnone WTnone_tpms.tsv
$ bash .command.sh 
[1] "Created by: riboviz"
[1] "Date: 2021-06-24 06:08:35"
...
  9.   +-vctrs:::stop_incompatible(...)
 10.     +-vctrs:::stop_vctrs(...)
Execution halted
```

**Note:** If rerunning a step that invokes a Python command then you need to provide a	`PYTHONPATH` with the path to your `riboviz` directory. For example:

```console
$ PYTHONPATH=$HOME/riboviz bash .command.sh
```

Alternatively, the command in `.command.sh` can be rerun directly, which may be desirable if you want to change input or configuration parameters when debugging, without rerunning the workflow via Nextflow. For example:

```console
$ Rscript --vanilla /home/ubuntu/riboviz/rscripts/collate_tpms.R \
  --tpms-file=TPMs_all_CDS_all_samples.tsv \
  WT3AT WT3AT_tpms.tsv WTnone WTnone_tpms.tsv
```

**Note:** If rerunning a step that invokes a Python command then you need to provide a	`PYTHONPATH` with the path to your `riboviz` directory. For example:

```console
$ PYTHONPATH=$HOME/riboviz python ...
```

For R you may want to run the script interactively, as described in the previous session, in which case you can run:

```console
$ Rscript --vanilla --args \
  --tpms-file=TPMs_all_CDS_all_samples.tsv \
  WT3AT WT3AT_tpms.tsv WTnone WTnone_tpms.tsv
```

then, within R, run:

```R
> source("../../../rscripts/collate_tpms.R")
```

providing the relative path to the R script, in this case, `collate_tpms.R`.

### Debugging within a directory emulating Nextflow's `work/` step-specific directories

Rather than using step-specific directories and their cryptic auto-generated name, you can mimic the way in which Nextflow runs steps by creating your own directory and populating it with the required input files (or symbolic links to these). For example, the following assumes that the workflow has run and one wishes to create a directory to debug `generate_stats_figs.R` using the `WTnone" sample:

* Create directory:

```console
$ mkdir debug_gen_stats_figs
$ cd debug_gen_stats_figs
```

* Create symbolic links to input files required by `generate_stats_figs.R` (these commands implicitly create symbolic links with the same names as the file names being linked to):

```console
$ ln -s ../vignette/output/WTnone/WTnone.h5
$ ln -s ../data/yeast_codon_pos_i200.RData
$ ln -s ../data/yeast_features.tsv 
$ ln -s ../data/yeast_standard_asite_disp_length.txt
$ ln -s ../data/yeast_tRNAs.tsv
$ ln -s ../vignette/input/yeast_YAL_CDS_w_250utrs.fa
$ ln -s ../vignette/input/yeast_YAL_CDS_w_250utrs.gff3
$ ls -1AF
WTnone.h5@
yeast_codon_pos_i200.RData@
yeast_features.tsv@
yeast_standard_asite_disp_length.txt@
yeast_tRNAs.tsv@
yeast_YAL_CDS_w_250utrs.fa@
yeast_YAL_CDS_w_250utrs.gff3@
```

* Run the command to be debugged within this directory. For example:

```console
$ Rscript --vanilla /home/ubuntu/riboviz/rscripts/generate_stats_figs.R ...
```

* Or, for an interactive session within R:

```console
$ R --vanilla --args ...
```
```R
> source("../rscripts/generate_stats_figs.R")
```
