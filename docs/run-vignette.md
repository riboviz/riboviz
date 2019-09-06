# Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file

RiboViz is designed to be run from a single script using a single configuration file in [YAML](http://www.yaml.org/) format containing all the information required for the analysis workflow.

`riboviz/tools/prep_riboviz.py` contains a Python implementation of the RiboViz analysis workflow.

Using this script, you can run a "vignette" of the workflow on a sample data set, to see RiboViz's capabilities. An example for *Saccharomyces cerevisiae* reads, up to the output of HDF5 files, is provided.

`vignette/vignette_config.yaml` contains the workflow's configuration information, in YAML format.

---

## Inputs

### Organism-specific data

`data/` holds genome and annotation data for *Saccharomyces cerevisiae* (yeast).

Sequence and annotation files for the whole approximate *Saccharomyces cerevisiae* transcriptome were prepared using [script_for_transcript_annotation.Rmd](../rmarkdown/script_for_transcript_annotation.Rmd). 

This script processes the files `saccharomyces_cerevisiae_R64-2-1_20150113.gff` and `S288C_reference_sequence_R64-2-1_20150113.fsa` from the file [genome release R64-2-1](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz) downloaded from the [Saccharomyces Genome Database](https://www.yeastgenome.org/). The script was used to produce the following files holding S288c annotations and ORF sequences:

* `data/yeast_CDS_w_250utrs.fa`: transcript sequences.
* `data/yeast_CDS_w_250utrs.gff3`: coding sequences locations.
* `data/yeast_codon_pos_i200.RData`: position of codons within each gene (the numbering ignores the first 200 codons).

`yeast_CDS_w_250utrs.fa` and `yeast_CDS_w_250utrs.gff3` were downsampled to provide a manageable data set for demonstration purposes. See [User data for analyses within the context of an organism](#user-data-for-analyses-within-the-context-of-an-organism) below.

When running RiboViz, `yeast_codon_pos_i200.RData` is read by [generate_stats_figs.R](../rscripts/generate_stats_figs.R) to help with generating plots and tables of results data.

### User data for analyses within the context of an organism

For an organism, the following data is required by RiboViz:

* Transcript sequences in a FASTA file. The transcript sequences need to contain both coding regions and flanking regions (which could be fixed, or coincident with measured UTRs).
* Locations of coding sequences within the transcripts in a GTF/GFF3 file.
* Ribosomal rRNA and other contaminant sequences in a FASTA file.

`vignette/input` holds sample data for *Saccharomyces cerevisiae*-related analyses.

**Downsampled yeast transcriptome**

As the yeast data files are very large, these were downsampled. The data files `yeast_CDS_w_250utrs.fa` and `yeast_CDS_w_250utrs.gff3` were processed by filtering only ORFs in the left arm of chromosome 1, for which the ORF names start with `YALnnnx`. This produced the yeast genome and annotation data files:

* `vignette/input/yeast_YAL_CDS_w_250utrs.fa`: transcript sequences to align to, from just the left arm of chromosome 1.
* `vignette/input/yeast_YAL_CDS_w_250utrs.gff3`: matched genome feature file, specifying coding sequences locations (start and stop coordinates).

`yeast_YAL_CDS_w_250utrs.fa` is an input specified within the `vignette/vignette_config.yaml` file's `orf_fasta` parameter.

`yeast_YAL_CDS_w_250utrs.gff3` is an input specified within the `vignette/vignette_config.yaml` file's `orf_gff_file` parameter.

The document [Appendix A1: Yeast Nomenclature Systematic Open Reading Frame (ORF) and Other Genetic Designations](https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527636778.app1) describes the ORF naming convention.

**Ribosomal RNA (rRNA) contaminants to remove**

Selected `RDN-n-n` sequences were copied and pasted from the file `rna_coding_R64-1-1_20110203.fasta` found in the file [genome release R64-1-1](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz) from the [Saccharomyces Genome Database](https://www.yeastgenome.org/) (note that this is one release earlier from that cited in [Organism-specific data](#organism-specific-data) above). The resulting file was:

* `vignette/input/yeast_rRNA_R64-1-1.fa`: rRNA sequences to avoid aligning to.

This file is an input specified within the `vignette/vignette_config.yaml` file's `rRNA_fasta` parameter.

**Vignette read data**

The read data files are downsampled ribosome profiling data from *Saccharomyces cerevisiae*. The data was downsampled to provide a dataset that was realistic, but small enough to run quickly.

The data is from the paper Guydosh N.R. and Green R. "[Dom34 rescues ribosomes in 3' untranslated regions](https://www.ncbi.nlm.nih.gov/pubmed/24581494)", Cell. 2014 Feb 27;156(5):950-62. doi: [10.1016/j.cell.2014.02.006](https://doi.org/10.1016/j.cell.2014.02.006). The NCBI accession for the whole dataset is #GSE52968.

* SRX386986: GSM1279570: wild-type no additive, [SRR1042855](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1042855).
* SRX386995: GSM1279579: wild-type plus 3-AT, [SRR1042864](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1042864).

In July 2017, these files were imported using NCBI's [fastq-dump](https://ncbi.github.io/sra-tools/fastq-dump.html) and gzipped to produce:

* `SRR1042855.fastq.gz`
* `SRR1042864.fastq.gz`

The data was sampled uniformly at random 1/50 reads from each file, producing about 1 million reads total, to produce the read data files:

* `vignette/input/SRR1042855_s1mi.fastq.gz`: ~1mi-sampled RPFs wild-type no additive.
* `vignette/input/SRR1042864_s1mi.fastq.gz`: ~1mi-sampled RPFs wild-type + 3-AT.

The data files are inputs specified within the `vignette/vignette_config.yaml` file's `fq_files` parameter:

```yaml
fq_files: # fastq files to be processed
  WTnone: SRR1042855_s1mi.fastq.gz # do not use "_" in dataset name
  WT3AT: SRR1042864_s1mi.fastq.gz
  NotHere: example_missing_file.fastq.gz # prep_riboviz should give error message for missing files
```

Note that the file specifies an additional, non-existent, file. This is used to test that the workflow processes valid files and ignores non-existent ones.

---

## What `prep_riboviz.py` does

The script prepares ribosome profiling data for RiboViz or other analyses. It does the following (`names` in brackets correspond to variables in the YAML configuration file):

* Reads configuration information from the YAML configuration file.
* Builds hisat2 indices if requested (`build_indices=TRUE`, by default), in an index directory (`dir_index`).
* Processes all `fastq.gz` files (`fq_files`) within the input directory  (`dir_in`).
* Cuts out sequencing library adapters (`adapters=CTGTAGGCACC`, by default) with `cutadapt`.
* Removes rRNA or other contaminating reads by hisat2 alignment to the rRNA index file (`rRNA_index`).
* Aligns remaining reads to ORFs or another hisat2 index file (`orf_index`).
* Trims 5' mismatches from reads and removes reads with more than 2 mismatches via a call to `riboviz/tools/trim5pmismatch.py`.
* Parallelizes over many processes (`nprocesses`):
  - This value is used to configure hisat2, samtools sort, bam_to_h5.R and generate_stats_figs.R.
  - For cutadapt and Python 3, the number of available processors on the host will be used.
  - For cutadapt and Python 2, its default of 1 processor will be used as cutadapt cannot run in parallel under Python 2.
* Makes length-sensitive alignments in compressed h5 format by running `rscripts/bam_to_h5.R`.
* Generates summary statistics, and analyses and QC plots for both RPF and mRNA datasets, by running `rscripts/generate_stats_figs.R`.
* Estimates read counts, reads per base, and transcripts per million for each ORF in each sample.
* Puts all intermediate files into a temporary directory (`dir_tmp`) which will be **large**.
* When finished, puts useful output files into an output directory (`dir_out`).
* Optionally exports bedgraph files for plus and minus strands (`make_bedgraph=TRUE`, by default).

For each sample or condition you want to compare (which should be placed into a single `.fastq.gz` file in the input directory (`dir_in`)), the configuration file needs a sub-variable of `fq_files`, whose name will be used in the output files. For example:

```yaml
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

---

## Warning: a temporary directory is created which could be very large!

Riboviz generates many intermediate files that could be large, i.e. about the same size as your input fasta files. All these files are placed in a temporary directory defined by a `dir_tmp` property in `vignette_config.yaml`, which, by default, is `vignette/tmp/`. Its contents can be inspected for troubleshooting if necessary. You should probably delete these temporary directories when you have completed your analysis.

For the vignette the total size of the temporary files is ~1141 MB, and the total size of all the files in the `vignette/` post-run is ~1152 MB.

---

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

* **Note:** `cutadapt` and `samtools`, which are invoked during the run, can only run under 1 process with Python 2.

### Run `prep_riboviz.py`

If you have not already done so, activate your Python environment:

* Miniconda Python 2.7+:

```console
$ source $HOME/miniconda2/bin/activate
```

* Miniconda Python 3.6+:

```console
$ source $HOME/miniconda3/bin/activate
```

If you have not already done so, set the paths to Hisat2 and Bowtie:

* If you followed [Create `setenv.sh` to configure Hisat2 and Bowtie paths](./install.md#create-setenvsh-to-configure-hisat2-and-bowtie-paths), then run:

```console
$ source $HOME/setenv.sh
```

* Otherwise, run the following (your directory names may be different, depending on the versions of Hisat2 and Bowtie you have):

```console 
export PATH=~/hisat2-2.1.0:$PATH
export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
```

Run `prep_riboviz.py`:

* Python 3:

```console
$ python -u -m riboviz.tools.prep_riboviz rscripts/ data/ \
    vignette/vignette_config.yaml 2>&1 | tee vignette-out.txt
```

* Python 2 or 3:

```console
$ PYTHONPATH=. python -u riboviz/tools/prep_riboviz.py rscripts/ data/ \
    vignette/vignette_config.yaml 2>&1 | tee vignette-out.txt
```

The part of the command, `2>&1 | tee vignette-out.txt`, allows for the messages that `prep_riboviz.py` prints to be both printed onto the screen and also printed into the `vignette-out.txt` file so you can inspect it later. This can be useful if any problems arise.

The exit code can be checked by running

```console
$ echo ${PIPESTATUS[0]}
```

If all went well an exit code of 0 will be returned.

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

```console
$ free --mega
              total        used        free      shared  buff/cache   available
Mem:           2017         684        1028           2         303        1181
Swap:           969         619         350
```

Divide the free memory by the number of processes, `nprocesses` e.g. 1024/4 = 256 MB.

Edit `riboviz/tools/prep_riboviz.py` and change the lines:

```python
    cmd_sort = ["samtools", "sort", "-@", str(nprocesses),
                "-O", "bam", "-o", fn_out + ".bam", "-"]
```

to include the `samtools` flag `-m <MEMORY_DIV_PROCESSES>M` e.g.:

```python
    cmd_sort = ["samtools", "sort", "-@", str(nprocesses),
                "-m", "256M",
                "-O", "bam", "-o", fn_out + ".bam", "-"]
```

### Check the expected output files

You should expect to see the following files produced.

Index files in `vignette/index`. For example:

```
YAL_CDS_w_250.1.ht2 # hisat2 indices from yeast_YAL_CDS_w_250utrs.fa
YAL_CDS_w_250.2.ht2
YAL_CDS_w_250.3.ht2
YAL_CDS_w_250.4.ht2
YAL_CDS_w_250.5.ht2
YAL_CDS_w_250.6.ht2
YAL_CDS_w_250.7.ht2
YAL_CDS_w_250.8.ht2
yeast_rRNA.1.ht2    # hisat2 indices from yeast_rRNA_R64-1-1.fa
yeast_rRNA.2.ht2
yeast_rRNA.3.ht2
yeast_rRNA.4.ht2
yeast_rRNA.5.ht2
yeast_rRNA.6.ht2
yeast_rRNA.7.ht2
yeast_rRNA.8.ht2
```

For this example, the index files occupy ~9 MB:

```console
$ du -sm vignette/index/
9	vignette/index/
```

Intermediate outputs in `vignette/tmp`. For example:

```
WT3AT_nonrRNA.fq         # trimmed non-rRNA reads
WT3AT_orf_map_clean.sam  # orf-mapped reads with mismatched nt trimmed
WT3AT_orf_map.sam        # orf-mapped reads
WT3AT_rRNA_map.sam       # rRNA-mapped reads
WT3AT_trim.fq            # trimmed reads
WT3AT_unaligned.sam      # unaligned reads

WTnone_nonrRNA.fq
WTnone_orf_map_clean.sam
WTnone_orf_map.sam
WTnone_rRNA_map.sam
WTnone_trim.fq
WTnone_unaligned.sam
```

**Note:** these are uncompressed and large. For this example, the intermediate files will occupy ~1040 MB:

```console
$ du -sm vignette/tmp/
1040	vignette/tmp/
```

Outputs in `vignette/output`. For example:

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

For this example, the output files occupy ~3 MB:

```console
$ du -sm vignette/output/
3	vignette/output/
```

---

## Cleaning up to run again

Before rerunning the vignette, delete the auto-generated index, temporary and output directories:

```console
$ rm -rf vignette/index
$ rm -rf vignette/tmp
$ rm -rf vignette/output
```

You might also want to do this if you have run the vignette with a missing R package, and then want to run it again from scratch. Alternatively, you might have edited the vignette and committed your changes, and be submitting a pull request.

---

## Using generated hisat2 indices

If you have already generated hisat2 indices for the same organism and annotation, set `build_indices` in the configuration file to `FALSE` and use the same index directory (`dir_index`) that you built in to.

---

## Customising the vignette

We suggest copying `vignette/vignette_config.yaml` and the rest of the `vignette` directory, and then customising it to fit your own datasets.

---

## Anatomy of `prep_riboviz.py`

A summary of the commands run and files output at each stage of a run of the vignette with the default configuration (but with `nprocesses: 4`) and input files.

### Inputs

Configuration file:

```
vignette/vignette_config.yaml
```

User input files in `vignette/input/`:

```
SRR1042855_s1mi.fastq.gz
SRR1042864_s1mi.fastq.gz
yeast_rRNA_R64-1-1.fa
yeast_YAL_CDS_w_250utrs.fa
yeast_YAL_CDS_w_250utrs.gff3
```

Organism data files in `data/`:

```
yeast_codon_pos_i200.RData
yeast_features.tsv
yeast_tRNAs.tsv
```

### Build indices for alignment

Build rRNA index: 

```
hisat2-build vignette/input/yeast_rRNA_R64-1-1.fa \
    vignette/index/yeast_rRNA
```

Outputs files to `vignette/index/`:

```
yeast_rRNA.1.ht2
yeast_rRNA.2.ht2
yeast_rRNA.3.ht2
yeast_rRNA.4.ht2
yeast_rRNA.5.ht2
yeast_rRNA.6.ht2
yeast_rRNA.7.ht2
yeast_rRNA.8.ht2
```

Build ORF index:

```
hisat2-build vignette/input/yeast_YAL_CDS_w_250utrs.fa \
    vignette/index/YAL_CDS_w_250
```

Outputs files to `vignette/index/`:

```
YAL_CDS_w_250.1.ht2
YAL_CDS_w_250.2.ht2
YAL_CDS_w_250.3.ht2
YAL_CDS_w_250.4.ht2
YAL_CDS_w_250.5.ht2
YAL_CDS_w_250.6.ht2
YAL_CDS_w_250.7.ht2
YAL_CDS_w_250.8.ht2
```

### Process file WT3AT (`SRR1042864_s1mi.fastq.gz`)

Cut illumina adapters i.e. cut out sequencing library adapters (CTGTAGGCACC or adapters):

```
cutadapt --trim-n -O 1 -m 5 -a CTGTAGGCACC -o vignette/tmp/WT3AT_trim.fq \
    vignette/input/SRR1042864_s1mi.fastq.gz -j 4
```

Outputs files to `vignette/tmp/`:

```
WT3AT_trim.fq
```

Map reads to RNA i.e. remove rRNA or other contaminating reads by HISAT2 alignment to rRNA index file:

```
hisat2 -p 4 -N 1 --un vignette/tmp/WT3AT_nonrRNA.fq \
    -x vignette/index/yeast_rRNA -S vignette/tmp/WT3AT_rRNA_map.sam \
    -U vignette/tmp/WT3AT_trim.fq
```

Outputs files to `vignette/tmp/`:

```
WT3AT_nonrRNA.fq
WT3AT_rRNA_map.sam
```

Map to ORFs with (mostly) default settings, up to 2 alignments i.e. align remaining reads to ORFs or other HISAT2 index file:

```
hisat2 -p 4 -k 2 --no-spliced-alignment --rna-strandness F --no-unal \
    --un vignette/tmp/WT3AT_unaligned.fq -x vignette/index/YAL_CDS_w_250 \
    -S vignette/tmp/WT3AT_orf_map.sam -U vignette/tmp/WT3AT_nonrRNA.fq
```

Outputs files to `vignette/tmp/`:

```
WT3AT_orf_map.sam
WT3AT_unaligned.fq
```

Trim 5' mismatched nt from reads and remove reads with more than 2 mismatches:

```
python riboviz/tools/trim_5p_mismatch.py -mm 2 \
    -in vignette/tmp/WT3AT_orf_map.sam \
    -out vignette/tmp/WT3AT_orf_map_clean.sam
```

Outputs files to `vignette/tmp/`:

```
WT3AT_orf_map_clean.sam
```

Convert sam (text) file to bam (compressed binary) file:

```
samtools view -b vignette/tmp/WT3AT_orf_map_clean.sam
```

Capture output piped from the above and sort bam file on genome:

```
samtools sort -@ 4 -O bam -o vignette/output/WT3AT.bam -
```

Outputs files to `vignette/output/`:

```
WT3AT.bam
```

Index bam file to create bai file:

```
samtools index vignette/output/WT3AT.bam
```

Outputs files to `vignette/output/`:

```
WT3AT.bam.bai
```

Calculate transcriptome coverage for plus strand and export as a bedgraph:

```
bedtools genomecov -ibam vignette/output/WT3AT.bam -bga -5 -strand +
```

streaming output into output file.

Outputs files to `vignette/output/`:

```
WT3AT_plus.bedgraph
```

Calculate transcriptome coverage for minus strand and export as a bedgraph:

```
bedtools genomecov -ibam vignette/output/WT3AT.bam -bga -5 -strand -
```

streaming output into output file.

Outputs files to `vignette/output/`:

```
WT3AT_minus.bedgraph
```

Make length-sensitive alignments in compressed h5 format:

```
Rscript --vanilla rscripts/bam_to_h5.R --Ncores=4 --MinReadLen=10 \
    --MaxReadLen=50 --Buffer=250 --PrimaryID=Name --SecondID=NULL \
    --dataset=vignette --bamFile=vignette/output/WT3AT.bam \
    --hdFile=vignette/output/WT3AT.h5 \
    --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 \
    --ribovizGFF=True --StopInCDS=False
```

Outputs files to `vignette/output/`:

```
WT3AT.h5
```

Generate summary statistics, analyses plots and QC plots for both RPF and mRNA datasets:

```
Rscript --vanilla rscripts/generate_stats_figs.R --Ncores=4 \
    --MinReadLen=10 --MaxReadLen=50 --Buffer=250 --PrimaryID=Name \
    --dataset=vignette --hdFile=vignette/output/WT3AT.h5 \
    --out_prefix=vignette/output/WT3AT \
    --orf_fasta=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True \
    --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 \
    --dir_out=vignette/output \
    --t_rna=data/yeast_tRNAs.tsv \
    --codon_pos=data/yeast_codon_pos_i200.RData \
    --features_file=data/yeast_features.tsv --do_pos_sp_nt_freq=True
```

This script uses `--dir_data=data/` to access and read `data/yeast_codon_pos_i200.RData` and `/data/yeast_tRNAs.tsv` also.

Outputs files to `vignette/output/`:

```
WT3AT_3nt_periodicity.pdf
WT3AT_3nt_periodicity.tsv
WT3AT_codon_ribodens.pdf
WT3AT_codon_ribodens.tsv
WT3AT_features.pdf
WT3AT_pos_sp_nt_freq.tsv
WT3AT_pos_sp_rpf_norm_reads.pdf
WT3AT_pos_sp_rpf_norm_reads.tsv
WT3AT_read_lengths.pdf
WT3AT_read_lengths.tsv
WT3AT_tpms.tsv
```

### Process file WTnone (`SRR1042855_s1mi.fastq.gz`)

The same commands are run but with `SRR1042855_s1mi.fastq.gz` being passed to `cutadapt`.

The same types of file are output, but with prefix `WTnone`.

### Collate TPMs across samples

Only successfully processed samples are collated.

```
Rscript --vanilla rscripts/collate_tpms.R \
    --dir_out=vignette/output \
    WTnone \
    WT3AT
```

Outputs files to `vignette/output/`:

```
TPMs_collated.tsv
```

### Exit codes

`prep_riboviz.py` returns the following exit codes:

* 0: Processing successfully completed.
* 1: Errors occurred loading configuration.
* 2: Error occurred during indexing.
* 3: No samples were provided.
* 4: No sample was processed successfully.
* 5: Error occurred during TPMs collation.
