# Debugging

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
> source('<RELATIVE_PATH_TO_RIBOVIZ_DIRECTORY>/rscripts/generate_stats_figs.R')
```

To debug a specific line of code, you could add a `browser()` statement in the source first. Alternatively, you could copy and paste the parts of the code you wanted to run, as long as earlier dependencies are run first (packages, importing command arguments, function definitions).

**Note:** at present, the RiboViz R scripts `bam_to_h5.R`, `generate_stats_figs.R` and `collate_tpms.R` import other RiboViz R scripts. If running these RiboViz R scripts interactively, via `R` and `source`, then the directory in which they are run must be such that `rscripts` is a sibling of an ancestor of the directory in which the script is run interactively. For example, running a script interactively within a sub-sub-directory of Nextflow's `work/` directory or a `debug_gen_stats_figs` directory (as described in the next section) where either of these directories are in the same directory as `rscripts`.

---

## Debugging and Nextflow

As described in [Nextflow `work/` directory](../user/prep-riboviz-operation.md#nextflow-work-directory) in [What the RiboViz workflow does]((../user/prep-riboviz-operation.md), Nextflow runs each task, or workflow step, in its isolated `work/` sub-sub-directory. All the input files the steps needs are provided in that directory (as symbolic links) and all the output files are written into that directory. Once the step is complete Nextflow copies (or "publishes") the created files into temporary or output directories as defined in the RiboViz configuration (`dir_tmp` and `dir_out` parameters). However, files in these temporary and output directories are not read by any of the workflow steps when the workflow is running. Think of of it as the steps having no knowledge of the temporary or output directories at all, only the files in their `work/` sub-sub-directory.

### Debugging within Nextflow's `work/` sub-sub-directories

These sub-sub-directories can be used for debugging steps and the Python and R scripts that implement these steps.

When Nextflow is run, each workflow step will be assigned a unique, auto-generated identifier (a so-called "hash"). For example:

```console
$ nextflow run prep_riboviz.nf -params-file \
    vignette/vignette_config.yaml -ansi-log false 
N E X T F L O W  ~  version 20.04.1
Launching `prep_riboviz.nf` [sharp_caravaggio] - revision: 62d3b58c5d
No such sample file (NotHere): example_missing_file.fastq.gz
[35/385067] Submitted process > buildIndicesrRNA (yeast_rRNA)
...
[b6/4d6e14] Submitted process > generateStatsFigs (WTnone)
[9a/aeb496] Submitted process > generateStatsFigs (WT3AT)
...
[ce/2368b2] Submitted process > countReads
Workflow finished! (OK)
```

In this example the hash for the step `generateStatsFigs (WTnone)` - the application of `generate_stats_figs.R` on the files for sample `WTnone` - is `b6/4d6e14`. This provides the location of the task's directory in the Nextflow `work/` directory:

```console
$ cd work/b6/4d6e1425e1f839cbc763b6d00f8a40/
```

This directroy includes the symbolic links to input files (denoted by `@`), output files, and Nextflow-generated command and log files:

```console
$ ls -1AF
3ntframe_bygene.tsv
3ntframe_propbygene.pdf
3nt_periodicity.pdf
3nt_periodicity.tsv
codon_ribodens.pdf
codon_ribodens.tsv
.command.begin
.command.err
.command.log
.command.out
.command.run
.command.sh
.exitcode
features.pdf
pos_sp_nt_freq.tsv
pos_sp_rpf_norm_reads.pdf
pos_sp_rpf_norm_reads.tsv
read_lengths.pdf
read_lengths.tsv
startcodon_ribogridbar.pdf
startcodon_ribogrid.pdf
tpms.tsv
WTnone.h5@
yeast_codon_pos_i200.RData@
yeast_features.tsv@
yeast_standard_asite_disp_length.txt@
yeast_tRNAs.tsv@
yeast_YAL_CDS_w_250utrs.fa@
yeast_YAL_CDS_w_250utrs.gff3@
```

The auto-generated file `.command.sh` has the exact command that Nextflow ran for the step:

```console
$ cat .command.sh 
#!/bin/bash -ue
Rscript --vanilla /home/ubuntu/riboviz/rscripts/generate_stats_figs.R
  --num-processes=1            --min-read-length=10
  --max-read-length=50            --buffer=250
  --primary-id=Name            --dataset=vignette
  --hd-file=WTnone.h5
  --orf-fasta-file=yeast_YAL_CDS_w_250utrs.fa            --rpf=true
  --output-dir=.            --do-pos-sp-nt-freq=true
  --t-rna-file=yeast_tRNAs.tsv
  --codon-positions-file=yeast_codon_pos_i200.RData
  --features-file=yeast_features.tsv
  --orf-gff-file=yeast_YAL_CDS_w_250utrs.gff3
  --asite-disp-length-file=yeast_standard_asite_disp_length.txt
  --count-threshold=64
```

The step can be rerun as follows:

```console
$ bash .command.sh
```

The outputs of the command will be written into the current directory.

Alternatively, the command in `.command.sh` can be rerun directly, which may be desirable if you want to change input or configuration parameters when debugging, without rerunning the workflow via Nextflow. For example:

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

For R you may want to run the script interactively, as described in the previous session, in which case you can run:

```console
$ Rscript --vanilla --args \
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

then, within R, run:

```R
> source("../../../rscripts/generate_stats_figs.R")
```

providing the relative path to the R script, in this case, `generate_stats_figs.R`.

### Debugging within a directory emulating Nextflow's `work/` sub-sub-directories

Rather than using `work/` sub-sub-directory and its cryptic auto-generated name, you can mimic the way in which Nextflow runs steps by creating your own directory and populating it with the required input files (or symbolic links to these). For example, the following assumes that the workflow has run and one wishes to create a directory to debug `generate_stats_figs.R` using the `WTnone" sample:

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
