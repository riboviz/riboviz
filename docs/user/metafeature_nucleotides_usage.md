# Running `metafeature_nucleotides.R`

`metafeature_nucleotides.R` is an R script that will produce a profile averaging the reads occurring around given nucleotide positions.

### Inputs 

`metafeature_nucleotides.R` requires 4 files as inputs:

* A H5 file - Produced by riboviz, lists all of the reads and their positions relative to the start codon.
* A GFF3 file - For an example species, this file can be found in `example-datasets` on the riboviz Github, in the `annotation` folder of the relevant species.
* A FASTA file - For an example species, this file can be found in `example-datasets` on the riboviz Github, in the `annotation` folder of the relevant species.
* A TSV file of nt positions to be studied - This file lists the gene under the heading 'Gene' and the nt postion under the heading 'Pos'. This is created by the user to fit their research question. For example: 

```
Gene	Pos
MIKE	5
MAT	7
MAT	5
MAT	6
```

### Outputs

`metafeature_nucleotides.R` will output a PDF containing a metafeature plot showing the average number of reads occurring at each position in a desired window around the position of interest. It is averaged over each position listed in the tsv file. 

### Arguments 

The arguments are: 


* '-i', '--input' - Path input to H5 file
* '-d', '--dataset' - Name of the dataset being studied
* '-g', '--gff' - Path to the GFF3 file of the organism being studied
* '-f', '--fasta' - Path to the FASTA file of the species being studied'),
* '--feature_pos' - A TSV file listing the Gene and Positions to normalize over, under the headings `Gene` and `Pos`
* '-o', '--output' - Path to output directory
* '--expand_width' - The desired range either side of the feature of interest, default = 5
* '--minreadlen' - The minimum read length, default = 10


### Command 

`Rscript rscripts/metafeature_nucleotides.R -i Mok-tinysim/output/A/A.h5 -d Mok-tinysim -g ../example-datasets/simulated/mok/annotation/tiny_2genes_20utrs.gff3 -f tiny_2genes_20utrs.fa --feature_pos [premade_tsv_file.tsv] -o .`
