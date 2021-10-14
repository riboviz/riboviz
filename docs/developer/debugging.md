# Debugging

* [Debugging R scripts with appropriate command-line arguments](#debugging-r-scripts-with-appropriate-command-line-arguments)
* [Debugging and Nextflow](#debugging and Nextflow)
  - [Debugging within Nextflow's `work/` step-specific directories](#debugging-within-nextflows-work-step-specific-directories)
  - [Debugging within a directory emulating Nextflow's `work/` step-specific directories](#debugging-within-a-directory-emulating-nextflows-work-step-specific-directories)

---

## Debugging R scripts with appropriate command-line arguments

To debug R scripts such as `generate_stats_figs.R` and `bam_to_h5.R`, they need to be run with the correct command-line arguments to discover the bug. R has good tools for interactive debugging, [explained in Hadley Wickham's chapter on debugging in R](https://adv-r.hadley.nz/debugging.html). However, interactive debugging tools such as `browser()` don't interrupt a call to `Rscript`. Instead you need to modify the call from:

```console
$ Rscript code_to_debug.R --myarg1 value1
```

to:

```console
$ R --args --myarg1 value1
```

then, from the R prompt run:

```R
> source('code_to_debug.R')
```

this will accept `debug()` and `browser()` statements run from the interactive R prompt.

For example, in the vignette we call:

```console
$ Rscript --vanilla /home/ubuntu/riboviz/rscripts/generate_stats_figs.R \
  --num-processes=1            --min-read-length=10 \
  --max-read-length=50            --buffer=250 \
  --primary-id=Name            --dataset=vignette \
  --hd-file=WTnone.h5 \
  --orf-fasta-file=yeast_YAL_CDS_w_250utrs.fa            --rpf=true \
  --output-dir=.            --do-pos-sp-nt-freq=true \
  --t-rna-file=yeast_tRNAs.tsv \
  --codon-positions-file=yeast_codon_pos_i200.RData \
  --features-file=yeast_features.tsv \
  --orf-gff-file=yeast_YAL_CDS_w_250utrs.gff3 \
  --asite-disp-length-file=yeast_standard_asite_disp_length.txt \
  --count-threshold=64
```

But to interactively debug a new feature, we'd run:

```console
$ R --vanilla --args \
  --num-processes=1            --min-read-length=10 \
  --max-read-length=50            --buffer=250 \
  --primary-id=Name            --dataset=vignette \
  --hd-file=WTnone.h5 \
  --orf-fasta-file=yeast_YAL_CDS_w_250utrs.fa            --rpf=true \
  --output-dir=.            --do-pos-sp-nt-freq=true \
  --t-rna-file=yeast_tRNAs.tsv \
  --codon-positions-file=yeast_codon_pos_i200.RData \
  --features-file=yeast_features.tsv \
  --orf-gff-file=yeast_YAL_CDS_w_250utrs.gff3 \
  --asite-disp-length-file=yeast_standard_asite_disp_length.txt \
  --count-threshold=64
```

then:

```R
> source('<PATH_TO_RIBOVIZ_DIRECTORY>/rscripts/generate_stats_figs.R')
```

To debug a specific line of code, you could add a `browser()` statement in the source first. Alternatively, you could copy and paste the parts of the code you wanted to run, as long as earlier dependencies are run first (packages, importing command arguments, function definitions).

**Note:** at present, the riboviz R scripts `bam_to_h5.R`, `generate_stats_figs.R` and `collate_tpms.R` import other riboviz R scripts. If running these riboviz R scripts interactively, via `R` and `source`, then the directory in which they are run must be such that `rscripts` is a sibling of an ancestor of the directory in which the script is run interactively. For example, running a script interactively within a sub-sub-directory of Nextflow's `work/` directory or a `debug_gen_stats_figs` directory (as described in the next section) where either of these directories are in the same directory as `rscripts`.

---

## Debugging and Nextflow

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
