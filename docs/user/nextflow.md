# RiboViz Nextflow workflow

[Nextflow](https://www.nextflow.io/) is a workflow management system. The file `prep_riboviz.nf` is a port of the RiboViz Python workflow, `riboviz.tools.prep_riboviz`, into a Nextflow workflow. In a future release, the RiboViz Python workflow will be deprecated by this Nextflow workflow.

---

## Configuring the Nextflow workflow

The Nextflow workflow uses the same YAML configuration file as `riboviz.tools.prep_riboviz` (as described in [Configuring the RiboViz workflow](./prep-riboviz-config.md).

---

## Running the Nextflow workflow

Run:

```console
$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE> -ansi-log false
```

where:

* `<CONFIG_FILE>`: path to a YAML configuration file.
* `-ansi-log false`: requests that each invocation of a Nextflow task is displayed on a separate line.

Configuration parameters can also be provided via the command-line in the form `--<PARAMETER>=<VALUE>` (for example `--make_bedgraph=FALSE`).

For example, to run the vignette:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [goofy_jones] - revision: 0e94b1fa62
No such file (NotHere): example_missing_file.fastq.gz
[93/e96307] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[f1/38397f] Submitted process > cutAdapters (WT3AT)
[1a/a01841] Submitted process > buildIndicesrRNA (yeast_rRNA)
[f0/5779d2] Submitted process > cutAdapters (WTnone)
[c2/497da3] Submitted process > hisat2rRNA (WTnone)
[14/479d3b] Submitted process > hisat2rRNA (WT3AT)
[cd/2e10ec] Submitted process > hisat2ORF (WTnone)
[02/1fef3f] Submitted process > trim5pMismatches (WTnone)
[9e/90df9d] Submitted process > samViewSort (WTnone)
[1e/1a3f4d] Submitted process > outputBams (WTnone)
[d8/3d85de] Submitted process > makeBedgraphs (WTnone)
[26/055ef0] Submitted process > bamToH5 (WTnone)
[d6/355a50] Submitted process > hisat2ORF (WT3AT)
[90/78ca31] Submitted process > trim5pMismatches (WT3AT)
[03/dfb679] Submitted process > samViewSort (WT3AT)
[85/ccbb35] Submitted process > outputBams (WT3AT)
[75/7e865e] Submitted process > bamToH5 (WT3AT)
[32/a55251] Submitted process > makeBedgraphs (WT3AT)
[d5/7a7e8b] Submitted process > generateStatsFigs (WTnone)
[f4/5188e3] Submitted process > generateStatsFigs (WT3AT)
Finished processing sample: WTnone
[c1/184176] Submitted process > renameTpms (WTnone)
Finished processing sample: WT3AT
[6f/716c9e] Submitted process > renameTpms (WT3AT)
[c6/bd5296] Submitted process > collateTpms (WTnone, WT3AT)
[a6/3cb086] Submitted process > countReads
Workflow finished! (OK)
```

Each sample-specific process is labelled with the sample ID, indexing processes are labelled with the index prefix, and multiplexed file-specific processes (not shown above, but an example is below) are labelled with the file name (minus extension)

The `collateTpms` process displays the names of all the samples that are collated.

---

## Help

Usage and configuration information can be viewed via use of the `--help` flag:

```console
$ nextflow run prep_riboviz.nf --help
```

Note that `--help` displays workflow help, whereas `-help` display's the `nextflow run` command's in-built help.

---

## Incremental build

If processing of a sample fails then `prep_riboviz.nf` has been written to ensure that this does not prevent the processing of other samples. For example, if the file for the vignette sample `WTnone` was corrupt in some way:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [mad_heisenberg] - revision: 134fcb8f6f
No such file (NotHere): example_missing_file.fastq.gz
[a4/e3fc5d] Submitted process > cutAdapters (WT3AT)
[03/f1cfc3] Submitted process > cutAdapters (WTnone)
[9d/43b3e0] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[ec/1aafe5] Submitted process > buildIndicesrRNA (yeast_rRNA)
[03/f1cfc3] NOTE: Process `cutAdapters (WTnone)` terminated with an error exit status (1) -- Error is ignored
[35/1c7ac8] Submitted process > hisat2rRNA (WT3AT)
[1a/af8ff1] Submitted process > hisat2ORF (WT3AT)
[c4/2f18a2] Submitted process > trim5pMismatches (WT3AT)
...
```

If a workflow fails, then a `-resume` option allows it to be rerun from the point at which it failed (for example, after fixing the issue with the file for `WTnone`):

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false -resume
```

This feature also supports incremental build. For example, given a `vignette_config.yaml` which specifies only sample `WTnone`, running Nextflow gives:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [serene_noether] - revision: 26bb98ec7b
No such file (NotHere): example_missing_file.fastq.gz
[bf/ddd7de] Submitted process > cutAdapters (WTnone)
[0d/f82601] Submitted process > buildIndicesrRNA (yeast_rRNA)
[61/5a1186] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[02/dac1f8] Submitted process > hisat2rRNA (WTnone)
[15/e5f9ed] Submitted process > hisat2ORF (WTnone)
[ed/fc50b9] Submitted process > trim5pMismatches (WTnone)
...
```

If `WT3AT` is then added to `vignette_config.yaml` and Nextflow is run with the `-resume` option, then only the processing for `WT3AT` is done, the cached outputs to date for `WTnone` being reused, not recomputed:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false -resume
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [desperate_coulomb] - revision: 26bb98ec7b
No such file (NotHere): example_missing_file.fastq.gz
[80/79a5f4] Submitted process > cutAdapters (WT3AT)
[bf/ddd7de] Cached process > cutAdapters (WTnone)
[61/5a1186] Cached process > buildIndicesORF (YAL_CDS_w_250)
[0d/f82601] Cached process > buildIndicesrRNA (yeast_rRNA)
[02/dac1f8] Cached process > hisat2rRNA (WTnone)
[15/e5f9ed] Cached process > hisat2ORF (WTnone)
[ed/fc50b9] Cached process > trim5pMismatches (WTnone)
[e3/9c39d7] Submitted process > hisat2rRNA (WT3AT)
[a6/49addd] Submitted process > hisat2ORF (WT3AT)
[46/d5db2c] Submitted process > trim5pMismatches (WT3AT)
...
```

---

## Multiplexed files

When running the workflow with a multiplexed FASTQ file (`multiplex_fq_files`) then a list of the sample IDs for which there were matching reads (and output files) and those for which there were no matching reads (and no output files) is printed.

For example, the sample sheet `data/simdata/multiplex_barcodes.tsv`, used in `vignette/simdata_multiplex_config.yaml` specifies a sample, `Tag3`, for which there will be no matching reads. The output of the workflow is as follows:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/simdata_multiplex_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [elated_jones] - revision: 0e94b1fa62
[e9/53999d] Submitted process > cutAdaptersMultiplex (multiplex_umi_barcode_adaptor)
[c3/6566db] Submitted process > buildIndicesrRNA (yeast_rRNA)
[43/752a9a] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[19/981102] Submitted process > extractUmisMultiplex (multiplex_umi_barcode_adaptor)
[e5/c52b82] Submitted process > demultiplex (multiplex_umi_barcode_adaptor)
Demultiplexed samples: [Tag0, Tag1, Tag2]
Non-demultiplexed samples: [Tag3]
[3f/73d0ab] Submitted process > hisat2rRNA (Tag2)
[18/c710a7] Submitted process > hisat2rRNA (Tag1)
[cc/f425d1] Submitted process > hisat2rRNA (Tag0)
[25/dd22c3] Submitted process > hisat2ORF (Tag2)
[b3/2db101] Submitted process > hisat2ORF (Tag1)
[85/9427ed] Submitted process > hisat2ORF (Tag0)
[f9/362265] Submitted process > trim5pMismatches (Tag0)
[d2/082ba5] Submitted process > trim5pMismatches (Tag2)
[de/e42cb0] Submitted process > trim5pMismatches (Tag1)
[63/1da506] Submitted process > samViewSort (Tag0)
[a5/ab590b] Submitted process > samViewSort (Tag2)
[72/4e9acd] Submitted process > samViewSort (Tag1)
[8b/c0936d] Submitted process > dedupUmis (Tag2)
[46/e9f43f] Submitted process > groupUmisPreDedup (Tag0)
[07/d696b3] Submitted process > groupUmisPreDedup (Tag2)
[52/2bdb04] Submitted process > dedupUmis (Tag0)
[36/ef62b5] Submitted process > dedupUmis (Tag1)
[5f/f51fca] Submitted process > groupUmisPreDedup (Tag1)
[9d/688af6] Submitted process > outputBams (Tag2)
[dc/e5d588] Submitted process > groupUmisPostDedup (Tag2)
[c5/690b4c] Submitted process > groupUmisPostDedup (Tag0)
[85/241fd3] Submitted process > outputBams (Tag0)
[9c/d88da8] Submitted process > bamToH5 (Tag2)
[d0/316e3d] Submitted process > makeBedgraphs (Tag2)
[a4/c2d620] Submitted process > makeBedgraphs (Tag0)
[0e/e327cb] Submitted process > bamToH5 (Tag0)
[62/789bb7] Submitted process > groupUmisPostDedup (Tag1)
[c7/13ea6b] Submitted process > outputBams (Tag1)
[8f/dabdca] Submitted process > bamToH5 (Tag1)
[d4/b4eef0] Submitted process > makeBedgraphs (Tag1)
[f8/6ba638] Submitted process > generateStatsFigs (Tag2)
[0a/06981b] Submitted process > generateStatsFigs (Tag1)
[f6/27751b] Submitted process > generateStatsFigs (Tag0)
Finished processing sample: Tag0
[52/447fbe] Submitted process > renameTpms (Tag0)
Finished processing sample: Tag1
[10/54496e] Submitted process > renameTpms (Tag1)
Finished processing sample: Tag2
[d7/7cc9ea] Submitted process > renameTpms (Tag2)
[6c/28e9e0] Submitted process > collateTpms (Tag0, Tag1, Tag2)
[29/354049] Submitted process > countReads
Workflow finished! (OK)
```

If no samples in the sample sheet have a matching read (and no sample-specific files are output at all) then the output will appear as follows:

```console
...
[53/c42278] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[47/59eaa1] Submitted process > buildIndicesrRNA (yeast_rRNA)
[58/24a15e] Submitted process > cutAdaptersMultiplex (multiplex_umi_barcode_adaptor)
[f8/df225f] Submitted process > extractUmisMultiplex (multiplex_umi_barcode_adaptor)
[53/8df442] Submitted process > demultiplex (multiplex_umi_barcode_adaptor)
Demultiplexed samples: []
Non-demultiplexed samples: [Tag0, Tag1, Tag2, Tag3]
Workflow finished! (OK)
```

---

## Debugging and bash scripts

As each `work/` subdirectory has symbolic links to any input files it requires,plus a bash script with the command that was run, that specific step can be run standalone. This can be useful for debugging. For example:

```console
$ cd work/37/b11ee1d2fb315a1b72adb65c11b44/
$ cat .command.sh 
#!/bin/bash -ue
hisat2 --version
hisat2 -p 1 -N 1 -k 1 --un nonrRNA.fq -x yeast_rRNA -S rRNA_map.sam -U trim.fq
$ bash .command.sh 
/home/ubuntu/hisat2-2.1.0/hisat2-align-s version 2.1.0
64-bit
Built on login-node03
Wed Jun  7 15:53:42 EDT 2017
Compiler: gcc version 4.8.2 (GCC) 
Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
952343 reads; of these:
  952343 (100.00%) were unpaired; of these:
    467194 (49.06%) aligned 0 times
    485149 (50.94%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
50.94% overall alignment rate
```

---

## Generating reports

`-with-report`, `-with-timeline` and `-with-dag` flags allow you to reques that Nextflow create reports on a run and an image of the task execution workflow. For example:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false \
    -with-report report.html -with-timeline timeline.html \
    -with-dag workflow.svg
```
