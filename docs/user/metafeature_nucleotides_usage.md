# Running `metafeature_nucleotides.R`

`metafeature_nucleotides.R` is an R script that will produce a profile averaging the reads occurring around given nucleotide positions.

## Inputs 

`metafeature_nucleotides.R` requires 4 files as inputs:

* A GFF3 file, for a species. Examples of these files can be found in the [example-datasets](https://github.com/riboviz/example-datasets) repository within the `annotation` folder of the relevant species.
* A FASTA file, for a species. Examples of these files can be found in the [example-datasets](https://github.com/riboviz/example-datasets) repository within the `annotation` folder of the relevant species.
* An H5 file, produced by riboviz, which lists all of the reads and their positions relative to the start codon. It is assumed that the H5 file was output during a riboviz run that used the above GFF3 and FASTA files.
* A TSV file of nucleotide (nt) positions to be studied. This file must list the gene under the heading 'Gene' and the nt postion under the heading 'Pos'. This is created by the user to fit their research question. For example: 

```
Gene	Pos
MIKE	5
MAT	7
MAT	5
MAT	6
```

## Outputs

`metafeature_nucleotides.R` will output a PDF containing a metafeature plot showing the average number of reads occurring at each position in a desired window around the position of interest. It is averaged over each position listed in the TSV file. 

## Usage

`metafeature_nucleotides.R` requires N arguments:

* `-i` or `--input`: Path to the H5 input file.
* `-d` or `--dataset`: Name of the dataset being studied.
* `-g` or `--gff`: Path to the GFF3 file of the organism being studied.
* `-f` or  `--fasta`: Path to the FASTA file of the species being studied.
* `--feature_pos`: A TSV file listing the Gene and Positions to normalize over, under the headings `Gene` and `Pos`.
* `-o` or `--output`: Path to output directory.

There are a number of optional arguments that can be used to change how `metafeature_nucleotides.R` processes the data. These are:

* `--expand_width`: The desired range either side of the feature of interest, default = 5
* `--minreadlen`: Minimum read length, default = 10.
* `--asite_length`: Path to species-specific A-site displacement length file, \
default = `data/yeast_standard_asite_disp_length.txt`.

`metafeature_nucleotides.R` can be run with the command:

```console
$ Rscript rscripts/metafeature_nucleotides.R
        -i [Path to H5 file] \
        -d [Dataset] \
        -g [Path to GFF3 file] \
        -f [Path to FASTA file] \
        --feature_pos [Path to feature positions file] \
        -o [Output directory]
```

## Examples

```console
$ Rscript rscripts/metafeature_nucleotides.R -i data/Mok-tinysim/A.h5 -d Mok-tinysim -g data/Mok-tinysim/tiny_2genes_20utrs.gff3 -f data/Mok-tinysim/tiny_2genes_20utrs.fa --feature_pos data/feature_pos.tsv -o . --minreadlen 10 --expand_width 2
```
