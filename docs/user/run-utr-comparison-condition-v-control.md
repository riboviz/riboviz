# Running utr_comparison_condition_v_control.R

`utr_comparison_condition_v_control.R` takes the reads from two samples and calculates the ratio of reads per base (rpb) in the 5'UTR compared to the rpb in the CDS for each transcript. 

## Input ##

The script takes two h5 files as inputs, which are compared against each other. Currently the script is desgined in a way where it is expected that both samples are from the same dataset. A GGF3 file is also necessary. 

## Output ##

The script outputs a scatter plot, with each point representing a different gene. The coordinates of the point are calculated by dividing the rpb of the 5UTR by the rpb of the CDS. The Y coordinate holds the 5'UTR:CDS use value for one sample and the X coordinate holds the 5'UTR:CDS use value for the other, typically the control. A point on or near the X=Y line does not have a large change in the 5UTR:CDS rpb value between samples.  

## Arguments ##

* "-i", "--input" - Path input to the first H5 file to be tested, which forms the Y axis  
* "-c", "--compare" - Path to the second H5 file, to be compared with input file. Ideally from a control
* "-d", "--dataset" - Name of the dataset being studied
* "-g", "--gff" - Path to the GFF3 file of the organism being studied
* "-o", "--output" - Path to output directory
* --gene" - Any genes of interest to be highlighted in plots (default is NULL)
* "-t", "--treatment" - Treatment used on the cells, for labelling of plots (default is NULL)
* "--read_threshold" - The minimum reads per base value needed in both samples for a gene to be included in the plot (default is 0.02 reads per base)

## Example ##

`Rscript rscripts/utr_comparison_condition_v_control.R -i B-Sc_2012/output/VEG_1/VEG_1.h5 -c B-Sc_2012/output/ANAPH_1/ANAPH_1.h5 -g ../example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cerevisiae_yeast_CDS_w_250utrs.gff3 -d B-Sc_2012`