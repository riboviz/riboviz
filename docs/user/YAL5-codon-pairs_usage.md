# Running YAL5-codon-pairs.R

## Input ##
The script `YAL5-codon-pairs.R` looks into metafeatures, currently the translational landscape surrounding codon pairs. 

- h5 file
- gff3 file
- tsv codon positions file

Read counts are fetched from the H5 file of a sample and assigns the counts to the cooresponding codon pair position for the genes using the species-specific GFF3 file and a tsv file containing the codon positions for all genes of the species being studied, i.e. the `yeast_codon_table.tsv`. The script calculates the normalised number of reads within a window surrounding each occurrence of the feature of interest (e.g. a codon pair) which are then overlayed to generate a single plot. Normalisation is carried out within a set window of codons adjacent to the feature of interest. The default value of this window is currently 5 codons either side of the feature of interest. 
The aim of this script is to study translation by examining whether specific features can slow down translation, thereby having an inhibitory effect. If the feature is being translated more slowly than its surrounding codons this will generate a peak. 

## Outputs ##
There is currently one output when running YAL5-codon-pairs.R and one that is being tested.

The first is a pdf of a graph plotting the average relative number of reads mapping to each relative position within the desired window of the feature of interest. This is produced when one feature is given, for example the codon pair 'CGA GCC'.

The second potential output that is currently being tested is a table in the form of a tsv file. This is produced when multiple features are given in the form of a tsv file; for example all possible codons. This table will contain the average relative count of the each feature, and rank them in descending order. This will allow for the identification of features that have a high average relative count, so are potentially inhibitory.

NOTE: When giving multiple features of interest, the script is expecting them to be in the first column of a tsv file.

## Arguments ##
Six arguments are required for this script (the rest are set as defaults):
- h5 file
- sample dataset
- gff3 file
- codon_table
- feature of interest 
- output directory 

The arguemnts that have set defaults can be changed to alter how the script processes the data. The size of the expanded window can be adjusted, reading frames can be filtered for (or not) but these are not required as they have default values.

As optparse hs not been incorporated into the YAL5-codon-pairs.R, the script can be run from the riboviz directory with this command from the console:
```
system('Rscript <PATH TO SCRIPT> -i <PATH TO H5 FILE> -d <DATASET> -g <PATH TO GFF3 FILE> -a <PATH TO CODON POSITIONS FILE> --feature <CODON PAIR> -o <PATH>')
```

Available arguments:
```
'-i', '--input'. Path input to h5 file
'-d', '--dataset'. Name of the dataset being studied
'-g', '--gff'. Path to the GFF3 file of the organism being studied
'-a', '--annotation'. Path to codon positions file for the organism (tsv file)
'--feature'. Feature of interest, e.g. codon pair 
'-o', '--output'. Path to output directory
'--expand_width'.The desired range either side of the feature of interest, default = 5L
'--startpos'. Position of the start codon, default = 1
'--startlen'. Smallest length of reads, default = 10
'--frame'. Logical - keep all reads or filter for reading frame, default = FALSE
'--minreadlen'. Minimum read length, default = 10
'--colsum_out'. Logical, default = TRUE
'--snapdisp'. Reading frame to filter for if `frame = TRUE`, default = 0L)
```
