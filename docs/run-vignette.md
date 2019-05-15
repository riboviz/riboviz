# Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file

RiboViz is designed to be run from a single script using a single configuration file in [YAML](http://www.yaml.org/) format containing all the information required for the analysis.

There are two implementations of this script available:

* `scripts/prepRiboviz.sh`, a bash shell implementation.
* `scripts/prepRiboviz.py`, a Python implementation.

Using either of these scripts, you can run a "vignette" of the back-end analysis pipeline, to demostrate RiboViz's capabilities. An example for *Saccharomyces cerevisiae* reads, up to the output of HDF5 files, is in the `vignette/` directory.

## `vignette` directory contents

`vignette/input/` contains input data:

```
yeast_YAL_CDS_w_250utrs.fa    # sequence file to align to from just left arm of chromosome 1.
yeast_YAL_CDS_w_250utrs.gff3  # matched genome feature file to specify start and stop co-ordinates.
yeast_rRNA_R64-1-1.fa         # rRNA sequence file to avoid aligning to.
SRR1042855_s1mi.fastq.gz      # about 1mi-sampled RPFs wild-type no additive from Guydosh & Green 2014.
SRR1042864_s1mi.fastq.gz      # about 1mi-sampled RPFs wild-type + 3-AT from Guydosh & Green 2014.
```

`vignette/vignette_config.yaml` contains configuration information in YAML. It specifies that the two `fastq.gz` files above are to be processed, along with an additional non-existent file which is used for testing:

```
fq_files: # fastq files to be processed
  WTnone: SRR1042855_s1mi.fastq.gz # do not use "_" in dataset name
  WT3AT: SRR1042864_s1mi.fastq.gz
  NotHere: example_missing_file.fastq.gz # prepRiboviz should give error message for missing files
```

## What `prepRiboviz.sh` and `prepRiboviz.py` do

Each script prepares ribosome profiling data for RiboViz or other analyses. They each do the following (`names` in brackets correspond to variables in the YAML configuration file):

* Reads configuration information from the YAML configuration file.
* Builds hisat2 indices if requested (`build_indices=TRUE`, by default), in an index directory (`dir_index`).
* Processes all `fastq.gz` files (`fq_files`) within the input directory  (`dir_in`).
* Cuts out sequencing library adapters (`adapters=CTGTAGGCACC`, by default) with `cutadapt`.
* Removes rRNA or other contaminating reads by hisat2 alignment to the rRNA index file (`rRNA_index`).
* Aligns remaining reads to ORFs or another hisat2 index file (`orf_index`).
* Trims 5' mismatches from reads and removes reads with more than 2 mismatches via a call to `scripts/trim5pmismatch.py`.
* Parallelizes over many processors (`nprocesses`), except for `cutadapt` which isn't parallel.
* Makes length-sensitive alignments in compressed h5 format by running `scripts/reads_to_list.R`.
* Generates summary statistics, and analyses and QC plots for both RPF and mRNA datasets, by running `scripts/generate_stats_figs.R`.
* Estimates read counts, reads per base, and transcripts per million for each ORF in each sample.
* Puts all intermediate files into a temporary directory (`dir_tmp`) which will be **large**.
* When finished, puts useful output files into an output directory (`dir_out`).
* Optionally exports bedgraph files for plus and minus strands (`make_bedgraph=TRUE`, by default).

For each sample or condition you want to compare (which should be placed into a single `.fastq.gz` file in the input directory (`dir_in`)), the configuration file needs a sub-variable of `fq_files`, whose name will be used in the output files. For example:

```
fq_files:
  WTnone: SRR1042855_s1mi.fastq.gz 
  WT3AT: SRR1042864_s1mi.fastq.gz
  Example: data.fastq.gz
```

For each of these names (e.g. `Example`), many output files are produced in the output directory (`dir_out`):

* `Example.bam`, bamfile of reads mapped to transcripts, can be directly used in genome browsers.
* `Example.bam.bai`, bam index file for `Example.bam`.
* `Example_minus.bedgraph`, (optional), bedgraph of reads from minus strand.
* `Example_plus.bedgraph`, (optional), bedgraph of reads from minus strand.
* `Example.h5`, length-sensitive alignments in compressed h5 format.
* `Example_3nt_periodicity.tsv`
* `Example_3nt_periodicity.pdf`
* `Example_read_lengths.tsv`
* `Example_read_lengths.pdf`
* `Example_pos_sp_nt_freq.tsv`
* `Example_pos_sp_rpf_norm_reads.pdf`
* `Example_pos_sp_rpf_norm_reads.tsv`
* `Example_features.pdf`
* `Example_tpms.tsv`
* `Example_codon_ribodens.tsv`
* `Example_codon_ribodens.pdf`

A summary file is also put in the output directory:

* `TPMs_collated.tsv`, tab-separated test file with the transcripts per million (tpm) for all samples 

For each of these names (e.g. `Example`), the intermediate files produced in the temporary directory (`dir_tmp`) are:

* `Example_trim.fq`, trimmed reads.
* `Example_nonrRNA.fq`, trimmed non-rRNA reads.
* `Example_rRNA_map.sam`, rRNA-mapped reads.
* `Example_orf_map.sam`, orf-mapped reads.
* `Example_orf_map_clean.sam`, orf-mapped reads with mismatched nt trimmed.
* `Example_unaligned.sam`, unaligned reads.

The `_unaligned.sam` files could be used to find common contaminants or translated sequences not in your orf annotation.

## Warning: a temporary directory is created which could be very large!

Riboviz generates many intermediate files that could be large, i.e. about the same size as your input fasta files. All these files are placed in a temporary directory defined by a `dir_tmp` property in `vignette_config.yaml`, which, by default, is `vignette/tmp/`. Its contents can be inspected for troubleshooting if necessary. You should probably delete these temporary directories when you have completed your analysis.

For the vignette the total size of the temporary files is ~1141 MB, and the total size of all the files in the `vignette/` post-run is ~1152 MB.

## Run the "vignette"

By default, RiboViz is configured to use 1 process. If using Python 3.x then you can configure RiboViz to use additional processes:

* Open `vignette/vignette_config.yaml` in an editor.
* Change:

```yaml
nprocesses: 1 # number of processes to parallelize over
```

* to:

```yaml
nprocesses: 4 # number of processes to parallelize over
```

* **Note:** `samtools`, which is invoked during the run, can only run under 1 process with Python 2.

### Run `scripts/prepRiboViz.sh`

Run:

```bash
bash scripts/prepRiboviz.sh vignette/vignette_config.yaml
```

### Run `scripts/prepRiboViz.py`

Run:

```bash
python scripts/prepRiboviz.py scripts/ vignette/vignette_config.yaml
```


### Troubleshooting: `File vignette/input/example_missing_file.fastq.gz not found`

If you see:

```
File vignette/input/example_missing_file.fastq.gz not found
```

then this is expected and can be ignored. The vignette includes an attempt to analyse a missing input file, for testing, which is expected to fail.

### Troubleshooting: `vignette/output/NotHere_tpms.tsv does not exist, returning empty list`

If you see:

```
In get_tpms(paste0(ddir, "/", fstem, fend), ORFs) :
  vignette/output/NotHere_tpms.tsv does not exist, returning empty list
```

then this is expected and can be ignored. This arises due to the missing input file described above which is used for testing.

### Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`

If using more than one process (`nprocesses` in `vignette/vignette_config.yaml` > 1) you might get the error:

```
samtools sort: couldn't allocate memory for bam_mem
```

This leads to a failure to create `.bam` and `.bam.bai` files are not in `vignette/output/`.

You may need to explicitly set the amount of memory per thread in calls to `samtools sort`.

Check how much free memory you have e.g.

```bash
free --mega
```
```
              total        used        free      shared  buff/cache   available
Mem:           2017         684        1028           2         303        1181
Swap:           969         619         350
```

Divide the free memory by the number of processes, `nprocesses` e.g. 1024/4 = 256 MB.

If using `prepRiboviz.sh` then edit `scripts/prepRiboviz.sh` and change the lines:

```bash
        echo samtools sort -@ ${nprocesses} -O bam -o ${fn_out}.bam -

        samtools sort -@ ${nprocesses} -O bam -o ${fn_out}.bam -
```

to include the `samtools` flag `-m <MEMORY_DIV_PROCESSES>M` e.g.:

```bash
        echo samtools sort -@ ${nprocesses} -m 256M -O bam -o ${fn_out}.bam -

        samtools sort -@ ${nprocesses} -m 256M -O bam -o ${fn_out}.bam -
```

If using `prepRiboviz.py` then edit `scripts/prepRiboviz.py` and change the lines:

```python
    cmd_sort = ["samtools", "sort", "-@", str(config["nprocesses"]),
                "-O", "bam", "-o", fn_out + ".bam", "-"]
```

to include the `samtools` flag `-m <MEMORY_DIV_PROCESSES>M` e.g.:

```python
    cmd_sort = ["samtools", "sort", "-@", str(config["nprocesses"]),
                "-m", "256M",
                "-O", "bam", "-o", fn_out + ".bam", "-"]
```

## Check the expected output files

You should expect to see the following files produced.

### Index files in `vignette/index`

```
YAL_CDS_w_250.*.ht2  # hisat2 index from yeast_YAL_CDS_w_250utrs.
yeast_rRNA.*.ht2     # hisat2 index from yeast_rRNA_R64-1-1.
```

For example:

```
YAL_CDS_w_250.1.ht2
YAL_CDS_w_250.2.ht2
YAL_CDS_w_250.3.ht2
YAL_CDS_w_250.4.ht2
YAL_CDS_w_250.5.ht2
YAL_CDS_w_250.6.ht2
YAL_CDS_w_250.7.ht2
YAL_CDS_w_250.8.ht2
yeast_rRNA.1.ht2
yeast_rRNA.2.ht2
yeast_rRNA.3.ht2
yeast_rRNA.4.ht2
yeast_rRNA.5.ht2
yeast_rRNA.6.ht2
yeast_rRNA.7.ht2
yeast_rRNA.8.ht2
```

For this example, the index files will occupy ~9 MB:

```bash
du -sm vignette/index/
```
```
9	vignette/index/
```

### Intermediate outputs in `vignette/tmp`

```
*_trim.fq            # trimmed reads.
*_nonrRNA.fq         # trimmed non-rRNA reads.
*_rRNA_map.sam       # rRNA-mapped reads.
*_orf_map.sam        # orf-mapped reads.
*_orf_map_clean.sam  # orf-mapped reads with mismatched nt trimmed.
*_unaligned.sam      # unaligned reads.
```

For example:

```
WT3AT_nonrRNA.fq
WT3AT_orf_map_clean.sam
WT3AT_orf_map.sam
WT3AT_rRNA_map.sam
WT3AT_trim.fq
WT3AT_unaligned.sam

WTnone_nonrRNA.fq
WTnone_orf_map_clean.sam
WTnone_orf_map.sam
WTnone_rRNA_map.sam
WTnone_trim.fq
WTnone_unaligned.sam
```

**Note:** these are uncompressed and large. For this example, the intermediate files will occupy ~1040 MB:

```bash
du -sm vignette/tmp/
```
```
1040	vignette/tmp/
```

### Outputs in `vignette/output`

For example:

```
TPMs_collated.tsv

WT3AT_3nt_periodicity.pdf
WT3AT_3nt_periodicity.tsv
WT3AT.bam
WT3AT.bam.bai
WT3AT_codon_ribodens.pdf
WT3AT_codon_ribodens.tsv
WT3AT_features.pdf
WT3AT.h5
WT3AT_minus.bedgraph
WT3AT_plus.bedgraph
WT3AT_pos_sp_nt_freq.tsv
WT3AT_pos_sp_rpf_norm_reads.pdf
WT3AT_pos_sp_rpf_norm_reads.tsv
WT3AT_read_lengths.pdf
WT3AT_read_lengths.tsv
WT3AT_tpms.tsv

WTnone_3nt_periodicity.pdf
WTnone_3nt_periodicity.tsv
WTnone.bam
WTnone.bam.bai
WTnone_codon_ribodens.pdf
WTnone_codon_ribodens.tsv
WTnone_features.pdf
WTnone.h5
WTnone_minus.bedgraph
WTnone_plus.bedgraph
WTnone_pos_sp_nt_freq.tsv
WTnone_pos_sp_rpf_norm_reads.pdf
WTnone_pos_sp_rpf_norm_reads.tsv
WTnone_read_lengths.pdf
WTnone_read_lengths.tsv
WTnone_tpms.tsv
```

For this example, the output files will occupy ~3 MB:

```bash
du -sm vignette/output/
```
```
3	vignette/output/
```

## Capturing output

To both display all output from the script that is printed at the terminal, and capture it into a file, run, for example:

```bash
bash scripts/prepRiboviz.sh vignette/vignette_config.yaml 2>&1 | tee vignette-script-bash.txt
```
```bash
python scripts/prepRiboviz.py scripts/ vignette/vignette_config.yaml  2>&1 | tee script-py.txt
```

## Cleaning up to run again

Before rerunning the vignette, delete the auto-generated index, temporary and output directories:

```bash
rm -rf vignette/index
rm -rf vignette/tmp
rm -rf vignette/output
```

You might also want to do this if you have run the vignette with a missing R package, and then want to run it again from scratch. Alternatively, you might have edited the vignette and committed your changes, and be submitting a pull request.

## Using generated hisat2 indices

If you have already generated hisat2 indices for the same organism and annotation, set `build_indices` in the configuration file to `FALSE` and use the same index directory (`dir_index`) that you built in to.

## Customising the vignette

We suggest copying `vignette/vignette_config.yaml` and the rest of the `vignette` directory, and then customising it to fit your own dataset.
