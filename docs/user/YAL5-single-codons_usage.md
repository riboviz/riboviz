# Running `YAL5-single-codons.R`

## Inputs 

`YAL5-single-codons.R` looks into metafeatures, currently specifically codon use. It takes reads from the H5 file of a sample and assigns them to the corresponding codon position in a gene using the species GFF3 file and a TSV file listing the codons for all genes of the species being studied, for example `data/yeast_codon_table.tsv`. `YAL5-single-codons.R` runs to calculate the average relative number of reads mapping to a feature of interest, relative to a window of codons adjacent to the feature. The default value of this window is 5 codons either side of the feature of interest. This can be used to study how inhibitory a feature is, as more inhibitory features will have a larger average relative read value. 

## Outputs

There are two potential outputs when running `YAL5-single-codons.R`.

The first potential output is a PDF of a graph plotting the average relative number of reads mapping to each position in the desired window around the feature of interest. This is produced if only one feature is givenas a feature of interest, for example codon 'CGA'.

The second potential output is a table in the form of a TSV file. This is produced when multiple features are given in the form of a TSV file as the input for the feature of interest argument; for example all possible codons. This table will contain the average relative count of each feature, and rank them in descending order. This will allow for the identification of features that have a high average relative count, so are potentially inhibitory. 

**Note:** When given a file with multiple features of interest, `YAL5-single-codons.R` expects them to be in the first column of the TSV file.

## Arguments

`YAL5-single-codons.R` requires six arguments:

* `-i` or `--input`: Path to the H5 input file.
* `-d` or `--dataset`: Name of the dataset of the sample being studied.
* `-g` or `--gff`: Path to the GFF3 file of the organism being studied.
* `-a` or `--annotation`: Path to codon table for the organism (TSV file).
* `--feature`: Feature of interest, e.g. 'CGA'. If multiple features (codons) are being studied then they should be contained within a TSV file.  
* `-o` or `--output`: Path to output directory.

There are a number of optional arguments that can be used to change how `YAL5-single-codons.R` processes the data, such as if a reading frame is filtered for and which one, but these are not required as they have default values:

* `--expand_width`: The desired range either side of the feature of interest to be used for normalization, default = 5.
* `--frame`:  Reading frame to be studied, default = 0.
* `--minreadlen`: Minimum read length, default = 10.
* `--filter_for_frame`: Counts for all reading frames per codon are summed and assigned to their corresponding codon. Keep all by not filtering (FALSE) or filter for the specific reading frame, specified in `--frame`, (TRUE), default = TRUE.
* `--snapdisp`: Reading frame to filter for, if `--filter_for_frame = TRUE`, default = 0.
* '--asite_length': Path to species-specific A-site displacement length file, default = `data/yeast_standard_asite_disp_length.txt`.

`YAL5-single-codons.R` can be run with the command:

```console
$ Rscript rscripts/YAL5-single-codons.R \
	-i [Path to H5 file] \
 	-d [Dataset] \
	-g [Path to GFF3 file] \
	-a [Path to codon table for organism] \
	 --feature [individual feature to be studied or path to TSV file containing multiple features] \
	-o [Output directory]
```

## Examples

The examples assume that `YAL5-single-codons.R` is being run from the riboviz folder.

### Single feature of interest

Running on data from the simulated dataset Mok-simYAL5 with the aim of investigating the codon pair 'CGA', the command would be:

```console
$ Rscript rscripts/YAL5-single-codons.R -i Mok-simYAL5/output/A/A.h5 -d Mok-simYAL5 -g ../example-datasets/simulated/mok/annotation/Scer_YAL_5genes_w_250utrs.gff3 -a data/yeast_codon_table.tsv --feature CGA -o .
```

Running `YAL5-single-codons.R` with a single feature of interest produces a PDF with the following image:

<img src="../images/Meta_feature_plot_CGA_Mok-simYAL5.JPG" alt="CGA Mok-simYAL5 meta feature plot" width="500"/>

### Multiple features of interest

Running on data from the simulated dataset Mok-simYAL5 with the aim of investigating multiple codons, provided as a TSV file, the command would be:

```console
$ Rscript rscripts/YAL5-single-codons.R -i Mok-simYAL5/output/A/A.h5 -d Mok-simYAL5 -g ../example-datasets/simulated/mok/annotation/Scer_YAL_5genes_w_250utrs.gff3 -a data/yeast_codon_table.tsv --feature data/codons.tsv -o .
```

Running `YAL5-single-codons.R` with a TSV file containing multiple features of interest produces a file containing the following output format:

```
Feature	RelCount
CGA	2.30033370411568
CTC	1.95209626966086
AGG	1.85391041724436
GCG	1.81696428571429
TGC	1.72677661461294
TGG	1.48950846086332
GCA	1.4872624105102
GTA	1.48627944639399
AGC	1.47948869857231
GTG	1.45266557210408
```
