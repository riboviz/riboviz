# Running `reads_on_individual_gene_transcript.R`

The script `reads_on_individual_gene_transcript.R` will produce a transcript profile for a gene of interest, mapping reads to positions along the transcript. Two profiles are produced, one showing total reads and the other showing reads per million. 

## Input ##

The files needed to run this script are an input h5 file, which is produced at as an output of the riboviz pipeline, and a gff3 file for the species being studied. The user will also need to know the gene name of the gene they want to study. Currwntly the script only takes one gene as an input.

## Output ##

`reads_on_individual_gene_transcript.R` produces two output pdfs. One showing the reads along the transcript in the form 'total reads' and the other in the form 'reads per million reads'.  

## Arguments ##

* '-i', '--input'- input H5 file, Path to the H5 file containing the data to be studied
* '-g', '--gff' - Path to GFF3 file corresponding to the species being studied
* '--gene' - The gene of interest, the name used should match the name listed in H5 and GFF3 files so it is a good idea to check
* '-d', '--dataset' - Name of the dataset being run, ie Mok-simYAL5. This must match the dataset listed in H5 file
* '-o', '--output' - Output directory for the PDFs

## Examples ##

`Rscript rscripts/reads_on_individual_gene_transcript.R -i Mok-simYAL5/output/A/A.h5 -g ../example-datasets/simulated/mok/annotation/Scer_YAL_5genes_w_250utrs.gff3 --gene YAL038W -d Mok-simYAL5 -o .`