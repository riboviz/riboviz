# Run a "vignette" of the back-end analysis pipeline

`vignette/` contains files to run a "vignette" of Riboviz up to the output of h5 files.

`vignette/input` contains input data:

```
yeast_YAL_CDS_w_250utrs.fa    # sequence file to align to from just left arm of chromosome 1.
yeast_YAL_CDS_w_250utrs.gff3  # matched genome feature file to specify start and stop co-ordinates.
yeast_rRNA_R64-1-1.fa         # rRNA sequence file to avoid aligning to.
SRR1042855_s1mi.fastq.gz      # about 1mi-sampled RPFs wild-type no additive from Guydosh & Green 2014.
SRR1042864_s1mi.fastq.gz      # about 1mi-sampled RPFs wild-type + 3-AT from Guydosh & Green 2014.
```

`vignette/vignette_config.yaml` contains configuration information.

## Run the "vignette"

If using Python 2.x then, edit `vignette/vignette_config.yaml`:

* Change:

```
nprocesses: 4 # number of processes to parallelize over
```

* to:

```
nprocesses: 1 # number of processes to parallelize over
```

* This needs to be done as `samtools`, which is invoked during the run, can only run under 1 process with Python 2.

Run:

```bash
bash scripts/prepRiboviz.sh vignette/vignette_config.yaml
```

You should expect to see the following files produced.

Index files, in `vignette/index`:

```
YAL_CDS_w_250.*.ht2  # hisat2 index from yeast_YAL_CDS_w_250utrs.
yeast_rRNA.*.ht2     # hisat2 index from yeast_rRNA_R64-1-1.
```

* For example:

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

Intermediate outputs, in `vignette/tmp`:

```
*_trim.fq            # trimmed reads.
*_nonrRNA.fq         # trimmed non-rRNA reads.
*_rRNA_map.sam       # rRNA-mapped reads.
*_orf_map.sam        # orf-mapped reads.
*_orf_map_clean.sam  # orf-mapped reads with mismatched nt trimmed.
*_unaligned.sam      # unaligned reads.
```

* **Note:** these are uncompressed and large.
* For example:

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

Outputs, in `vignette/output`:

* For example:

```
WTnone.bam
WTnone.bam.bai
WTnone_minus.bedgraph
WTnone_plus.bedgraph
WTnone.h5
WTnone_3nt_periodicity.tsv
WTnone_3nt_periodicity.pdf
WTnone_read_lengths.tsv
WTnone_read_lengths.pdf
WTnone_pos_sp_nt_freq.tsv
WTnone_pos_sp_rpf_norm_reads.pdf
WTnone_pos_sp_rpf_norm_reads.tsv
WTnone_features.pdf
WTnone_features.tsv
WTnone_codon_ribodens.tsv
WTnone_codon_ribodens.pdf

WT3AT.bam
WT3AT.bam.bai
WT3AT_minus.bedgraph
WT3AT_plus.bedgraph
WT3AT.h5
WT3AT_3nt_periodicity.tsv
WT3AT_3nt_periodicity.pdf
WT3AT_read_lengths.tsv
WT3AT_read_lengths.pdf
WT3AT_pos_sp_nt_freq.tsv
WT3AT_pos_sp_rpf_norm_reads.pdf
WT3AT_pos_sp_rpf_norm_reads.tsv
WT3AT_features.pdf
WT3AT_features.tsv
WT3AT_codon_ribodens.tsv
WT3AT_codon_ribodens.pdf
```

## Before re-running

Before rerunning the vignette, delete the auto-generated files and folders:

```
run -rf vignette/index
run -rf vignette/tmp
run -rf vignette/output
```

## Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`

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

Divide the free memory by the number of processes, `nprocesses` e.g. 1024/4 = 256MB.

Edit `scripts/prepRiboviz.sh` and change the lines:


```
        echo samtools sort -@ ${nprocesses} -O bam -o ${fn_out}.bam -

        samtools sort -@ ${nprocesses} -O bam -o ${fn_out}.bam -
```

to include the `samtools` flag `-m <MEMORY_DIV_PROCESSES>M` e.g.:

```
        echo samtools sort -@ ${nprocesses} -m 256M -O bam -o ${fn_out}.bam -

        samtools sort -@ ${nprocesses} -m 256M -O bam -o ${fn_out}.bam -
```

## Troubleshooting: capturing output

To both display all output from the script that is printed at the terminal, and capture it into a file, run:

```bash
bash scripts/prepRiboviz.sh vignette/vignette_config.yaml 2>&1 | tee vignette-script.txt
```
