# Running `metafeature_nucleotides.R`

`metafeature_nucleotides.R` is an R script that will produce a profile averaging the reads occuring around given nucleotide positions.

### Inputs 

`metafeature_nucleotides.R` requires 4 files as inputs:

* A h5 file - Produced by riboviz, lists all of the reads and their positions relative to the start codon.
* A gff file - For an example species, this file can be found in `example-datasets` on the riboviz github, in the `annotation` folder of the relevant species.
* A fasta file - For an example species, this file can be found in `example-datasets` on the riboviz github, in the `annotation` folder of the relevant species.
* A tsv file of nt positions to be studied - This file lists the gene under the heading 'Gene' and the nt postion under the heading 'Pos'. This is created by the user to fit their research question.

### Outputs

`metafeature_nucleotides.R` will output a PDF containing a metafeature plot showing the avergage number of reads occuring at each position in a desired window around the position of interest. It is averaged over each position listed in the tsv file. 

### Arguments 

The optional arguments are 

```
option_list <- list(make_option(c('-i', '--input'),type = "character", help='Path input to h5 file'),
                    make_option(c('-d', '--dataset'),type = "character", help='Name of the dataset being studied'),
                    make_option(c('-g', '--gff'),type = "character", help='Path to the GFF3 file of the organism being studied'),
                    make_option(c('-f', '--fasta'),type = "character", help='Path to the fasta file of the organsim being studied'),
                    make_option(c('--feature_pos'), type = "character", help='A TSV file listing the Gene and Positions to normalize over'),
                    make_option(c('-o', '--output'), type = "character", help='Path to output directory'),
                    make_option(c('--expand_width'), type = "integer", help='the desired range either side of the feature of interest', default = 5),
                    make_option(c('--minreadlen'),type = "integer", help='minimum read length', default = 10))

```

### Command 

`Rscript rscripts/metafeature_nucleotides.R -i Mok-tinysim/output/A/A.h5 -d Mok-tinysim -g ../example-datasets/simulated/mok/annotation/tiny_2genes_20utrs.gff3 -f tiny_2genes_20utrs.fa --feature_pos [premade_tsv_file.tsv] -o .`