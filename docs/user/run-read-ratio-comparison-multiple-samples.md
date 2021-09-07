### Run read_ratio_comparison_multiple_samples.R

`read_ratio_comparison_multiple_samples.R` is an Rscript which compares the variance of reads mapping to the UTR relative to those mapping to those in the CDS between multiple samples from the same dataset. This script requires the R package `argparse` as multiple samples are given as part of the same argument.

## inputs 

The required input files are:

* A h5 file for each samole being compared 
* A GFF3 for the species being studied 

## Outputs

The script outputs a PDF containing a graph which has a number of boxplots equal to the number of samples used as an input. Each box shows the spread of UTR read per base : CDS reads per base in that sample. Each blue dot represents a gene.

## Arguments

* '-i', '--input' - Path to the H5 file containing the data to be studied. Able to take multiple files, with file paths separated by a space
* '-g', '--gff' - Path to GFF3 file corresponding to the species being studied
* '--gene'- Gene of interest, name used should match the name listed in H5 file so check format, default = NULL
* '-d', '--dataset'- Name of the dataset being run, ie D-Sp_2018, should match the dataset listed in H5 file
* '-o', '--output' - Output directory for plots, default = "."
* '--read_threshold' - The minimum number of reads per base needed for a gene to be plotted, default=0.02

## Example 

`Rscript rscripts/read_ratio_comparison_multiple_samples.R -i B-Sc_2012/output/VEG_1/VEG_1.h5 B-Sc_2012/output/ANAPH_1/ANAPH_1.h5 -g ../example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cerevisiae_yeast_CDS_w_250utrs.gff3 -d B-Sc_2012 -o .`