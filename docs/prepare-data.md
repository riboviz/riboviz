# Prepare data

For your organism, you will need:

* Transcript sequences in a fasta file. The transcript sequences need to contain both coding regions and flanking regions (which could be fixed, or coincident with measured UTRs). 
* Ribosomal rRNA and other contaminant sequences in a fasta file. 
* Locations of coding sequences within the transcripts in a GTF/gff3 file.

## Example: download transcript and rRNA sequences and gff3 file

For example, for_S. cerevisiae_, download:

* Transcript sequences, [yeast_CDS_w_250utrs.fa](https://raw.githubusercontent.com/RiboViz/RiboViz/master/data/yeast_CDS_w_250utrs.fa) (downloads from this Git repository)
* rRNA sequences, [rrna.fa](http://gdurl.com/KGnn/download)
* Coding sequences locations, [yeast_CDS_w_250utrs.gff3](https://raw.githubusercontent.com/RiboViz/RiboViz/master/data/yeast_CDS_w_250utrs.gff3) (downloads from this Git repository)

To download these files via bash:

```bash
mkdir sample-data
wget https://raw.githubusercontent.com/RiboViz/RiboViz/master/data/yeast_CDS_w_250utrs.fa -P sample-data
wget http://gdurl.com/KGnn/download -O sample-data/rrna.fa
wget https://raw.githubusercontent.com/RiboViz/RiboViz/master/data/yeast_CDS_w_250utrs.gff3 -P sample-data
```

The file sizes, and space requirements, are as follows:

```
5.6K	rrna.fa
 12M	yeast_CDS_w_250utrs.fa
885K	yeast_CDS_w_250utrs.gff3
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

It is easiest if you put all these files in a single input folder. Alternatively, you could symlink them from a single folder. For example, given:

```
samples/
    condition1/
        condition1.fastq.gz
    condition2/
        condition2.fastq.gz
```

This can be symlinked into a single folder as follows:

```bash
mkdir data
cd data/
ln -s ../samples/condition1/condition1.fastq.gz
ln -s ../samples/condition2/condition2.fastq.gz
ls -l
```
```
total 0
lrwxrwxrwx 1 ubuntu ubuntu 41 Aug 22 02:55 condition1.fastq.gz -> ../samples/condition1/condition1.fastq.gz
lrwxrwxrwx 1 ubuntu ubuntu 41 Aug 22 02:55 condition2.fastq.gz -> ../samples/condition2/condition2.fastq.gz
```
