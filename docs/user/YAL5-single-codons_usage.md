# Running `YAL5-single-codons.R`

## Inputs 

`YAL5-single-codons.R` looks into metafeatures, currently specifically codon use. It takes reads from the H5 file of a sample and assigns them to the corresponding codon position in a gene using the species GFF3 file and a TSV file listing the codons for all genes of the species being studied, for example `data/yeast_codon_table.tsv`. `YAL5-single-codons.R` runs to calculate the average relative number of reads mapping to a feature of interest, relative to a window of codons adjacent to the feature. The default value of this window is 5 codons either side of the feature of interest. This can be used to study how inhibitory a feature is, as more inhibitory features will have a larger average relative read value. 

## Outputs

There are two potential outputs when running `YAL5-single-codons.R`.

The first potential output is a PDF of a graph plotting the average relative number of reads mapping to each position in the desired window around the feature of interest. This is produced if only one feature is givenas a feature of interest, for example codon 'CGA'.

The second potential output is a table in the form of a TSV file. This is produced when multiple features are given in the form of a TSV file as the input for the feature of interest argument; for example all possible codons. This table will contain the average relative count of each feature, and rank them in descending order. This will allow for the identification of features that have a high average relative count, so are potentially inhibitory. 

**Note:** When given a file with multiple features of interest, `YAL5-single-codons.R` expects them to be in the first column of the TSV file.

## Execution

`YAL5-single-codons.R` requires six arguments:

* `-i` or `--input`: Path to the H5 input file.
* `-d` or `--dataset`: Name of the dataset of the sample being studied.
* `-g` or `--gff`: Path to the GFF3 file of the organism being studied.
* `-a` or `--annotation`: Path to codon table for the organism (TSV file).
* `--feature`: Feature of interest, e.g. 'CGA'. If multiple features (codons) are being studied then they should be contained within a TSV file with a 'Codon' header.
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

```console
$ cat Feature_Relative_use_Mok-simYAL5.tsv 
Feature	RelCount
CGA	4.90341753343239
GCG	4.42843465773097
TGG	2.87182281220301
GGG	2.87141812865497
TGC	2.85005325770739
AGG	2.81812422335214
GCA	2.48130673065506
AGC	2.17250602239257
GTG	2.16112635435094
ACG	2.0910912790197
ACA	2.09051734274515
GGA	2.06182344941109
GAG	2.03833092628914
CAG	1.99130324025171
GGC	1.94372412313873
CCC	1.90158625020691
CCA	1.84426110176067
CTC	1.83470120013005
CTG	1.7582253792347
AGA	1.75157049837252
GGT	1.66427555443987
TCG	1.61199993237037
AAG	1.60325385402508
GAC	1.46432952772166
CAT	1.3905267327226
TAC	1.33035411741904
ATA	1.32962685636946
AAA	1.31941516925201
GAT	1.27536928155255
TCA	1.25690230265607
CAC	1.22494298045999
CCT	1.18697580137237
GAA	1.13613286821688
GTA	1.04067550798979
CTA	1.01142881086673
CAA	0.97459463847631
ATG	0.956673097463448
CGT	0.94852950776106
CTT	0.932199733193517
AGT	0.897286292205183
TAT	0.763544818583199
TCC	0.696949289512051
TGT	0.683190675993757
TTG	0.624428138275046
ACC	0.553290563707026
GCT	0.533870191441569
ATC	0.49743367188103
AAC	0.494619421473075
ACT	0.463666867085569
TTA	0.451089555381943
GCC	0.383408468414928
AAT	0.378307702575587
TTC	0.37650940273614
GTC	0.305048109753676
TTT	0.27776876778166
TCT	0.27413728182794
ATT	0.0331192331093087
GTT	0
```

**Note:** If the TSV file has only one codon then a PDF is output as for the Single feature of interest mode of operation.
