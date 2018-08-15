# Prepare data

For your organism, you will need:

* Transcript sequences in a fasta file. The transcript sequences need to contain both coding regions and flanking regions (which could be fixed, or coincident with measured UTRs). 
* Ribosomal rRNA and other contaminant sequences in a fasta file. 
* Locations of coding sequences within the transcripts in a GTF/gff3 file.

## Example: download transcript and rRNA sequences and gff3 file

For example, for_S. cerevisiae_, download:

* Transcript sequences, [yeast_CDS_w_250utrs.fa](https://github.com/shahpr/RiboViz/blob/master/scripts/yeast_CDS_w_250utrs.fa)
* rRNA sequences, [rrna.fa](http://gdurl.com/KGnn/download)
* Coding sequences locations, [yeast_CDS_w_250utrs.gff3](https://github.com/shahpr/RiboViz/blob/master/scripts/yeast_CDS_w_250utrs.gff3)

To download these files via bash:

```bash
mkdir sample-data
wget https://github.com/shahpr/RiboViz/blob/master/scripts/yeast_CDS_w_250utrs.fa -P sample-data
wget http://gdurl.com/KGnn/download -O sample-data/rrna.fa
wget https://github.com/shahpr/RiboViz/blob/master/scripts/yeast_CDS_w_250utrs.gff3 -P sample-data
```

## Prepare input files

For each condition/sample in your experiment, merge all fastq files into a single gzipped fastq file, e.g. `condition1.fastq.gz`. The files are gzip-compressed to save space. For example:

```bash
cat condition1_subsetA.fastq condition1_subsetB.fastq | gzip >> condition1.fastq.gz
cat condition2_subsetC.fastq condition2_subsetD.fastq | gzip >> condition2.fastq.gz
.
.
.
```

It is easiest if you put all these files in a single input folder. Alternatively, you could symlink them from a single folder.
