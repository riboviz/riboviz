# riboviz output files and figures

This document lists the output files from a typical riboviz run, along with short descriptions. When first looking at outputs of a riboviz run, `read_counts_by_length.pdf` and `metagene_start_stop_read_counts.pdf` will be a good starting point for accessing the success of the run and the quality of the data used, as described below. 

**Initial Quality control checks**

When looking at the `read_counts_by_length.pdf`, it is expected to see a peak of reads at 28-31nt. If the adapter or UMI (Unique Molecular Identifiers) sequence given in the config file does not match those used in the experiment irregularities will be present in the output files. If these sequences are incorrect a clear peak at 28-31nt in `read_counts_by_length.pdf` will not be present. In addition, very few reads may be present as they are unable to align to the annotation files so are discarded at an earlier stage of the pipeline. You can examine `read_counts_per_file.tsv` to check the number of aligned reads at each stage of the pipeline. 

When looking at the `metagene_start_stop_read_counts.pdf`, it is expected to see a peak of reads upstream of the start codon, then repeating smaller peaks every 3nt. Severe disruption of the 3nt periodicity indicates a lower quality dataset and potential issues with the adapter sequence or UMI regular expression. 


Descriptions of specific files can be located quickly using the following links:

**Output files for an entire run**
* [TPMs_all_CDS_all_samples](#tpms_all_cds_all_samplestsv)
* [read_counts_per_file.tsv](#read_counts_per_filetsv)

**Output files for each samples** 
* [<SAMPLE_ID>_output_report.html](#sample_id_output_reporthtml)
* [<SAMPLE_ID>.bam](#sample_idbam)
* [<SAMPLE_ID>.bam.bai](#sample_idbambai)
* [minus.bedgraph](#minusbedgraph)
* [plus.bedgraph](#plusbedgraph)
* [<SAMPLE_ID>.h5](#sample_idh5)
* [metagene_start_stop_read_counts.tsv](#metagene_start_stop_read_countstsv)
* [metagene_start_stop_read_counts.pdf](#metagene_start_stop_read_countspdf)
* [read_counts_by_length.tsv](#read_counts_by_lengthtsv)
* [read_counts_by_length.pdf](#read_counts_by_lengthpdf)
* [pos_sp_nt_freq.tsv](#pos_sp_nt_freqtsv)
* [metagene_normalized_profile_start_stop.pdf](#metagene_normalized_profile_start_stoppdf)
* [metagene_normalized_profile_start_stop.tsv](#metagene_normalized_profile_start_stoptsv)
* [ORF_TPMs_vs_features.tsv](#orf_tpms_vs_featurestsv)
* [ORF_TPMs_vs_features.pdf](#orf_tpms_vs_featurespdf)
* [tpms.tsv](#tpmstsv)
* [normalized_density_APEsites_per_codon.tsv](#normalized_density_APEsites_per_codontsv)
* [normalized_density_APEsites_per_codon_long.tsv](#normalized_density_APEsites_per_codon_longtsv)
* [normalized_density_APEsites_per_codon.pdf](#normalized_density_APEsites_per_codonpdf)
* [gene_position_length_counts_5start.tsv](#gene_position_length_counts_5starttsv)
* [metagene_start_barplot_by_length.pdf](#metagene_start_barplot_by_lengthpdf)
* [metagene_start_ribogrid_by_length.pdf](#metagene_start_ribogrid_by_lengthpdf)
* [read_frame_per_ORF.tsv](#read_frame_per_ORFtsv)
* [read_frame_per_ORF_filtered.tsv](#read_frame_per_ORF_filteredtsv)
* [frame_proportions_per_ORF.pdf](#frame_proportions_per_ORFpdf)

After a riboviz run, many output files are produced within the output directory.
The output directory is specified by the parameter `dir_out` in `config.yaml`.

There are a few output files that collect information for an entire run.
There are many output files that are specific to each sample, which are organized into a separate subdirectory for each sample. 


# Output files for an entire run


## `TPMs_all_CDS_all_samples.tsv` 

A tsv file with transcripts per million (tpm) for all genes from successfully processed samples. This file is produced by the script `collate_TPMs.R`. This script uses the `tpms.tsv` file from each processed sample and lists the tmps of each sample for each gene, allowing for comparison between samples. 

```
ORF	ANAPH_1	VEG_1
Q0045	54.5	0.7
Q0050	0.6	0.7
Q0055	0.6	0.1
Q0060	2.2	0.1
Q0065	22.2	1
Q0070	22.8	0.7
```


## `read_counts_per_file.tsv` 

A [read counts file](#read-counts-file) (only if `count_reads: TRUE`).

A tsv file produced by `count_reads.py`. This file lists the sample being processed, the stage of the riboviz process, the path to the file being processed at that stage, the number of reads present at each stage, and includes a description of the process. The number of reads is expected to decrease between stages as reads are filtered and trimmed.

```
SampleName	Program	File	NumReads	Description
VEG_1	input	/exports/csce/eddie/biology/groups/wallace_rna/riboviz-emma/riboviz/riboviz/B-Sc_2012/input/SRR387871.fastq.gz	21155927	input
ANAPH_1	input	/exports/csce/eddie/biology/groups/wallace_rna/riboviz-emma/riboviz/riboviz/B-Sc_2012/input/SRR387890.fastq.gz	10210064	input
ANAPH_1	cutadapt	/exports/csce/eddie/biology/groups/wallace_rna/riboviz-emma/riboviz/riboviz/B-Sc_2012/tmp/ANAPH_1/trim.fq	9458608	Reads after removal of sequencing library adapters
ANAPH_1	hisat2	/exports/csce/eddie/biology/groups/wallace_rna/riboviz-emma/riboviz/riboviz/B-Sc_2012/tmp/ANAPH_1/nonrRNA.fq	5999936	rRNA or other contaminating reads removed by alignment to rRNA index files
```


# Output files for each sample

For each sample (`<SAMPLE_ID>`), intermediate files are produced in a sample-specific subdirectory (`<SAMPLE_ID>`).


## `<SAMPLE_ID>_output_report.html` 

This output report in .html format contains a provenance section for the sample, all the pdf output files produced by riboviz, and information on any plots not produced such as which files would be needed to produce plots in future runs. The HTML also includes a side bar to allow for navigation between figures. This file is produced by `AnalysisOutputs.Rmd`, which loads and creates all of the output graphs in HTML format.

Example of the top of the HTML file:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/html_output.jpg" alt="html_output image" width="500"/>

Only output if `run_static_html: TRUE`. 


## `<SAMPLE_ID>.bam` 

BAM file of reads mapped to transcripts, which can be directly used in genome browsers. The BAM file is produced by samtools.


## `<SAMPLE_ID>.bam.bai` 

BAM index file for `<SAMPLE_ID>.bam`. This file is produced by samtools.


## `minus.bedgraph` 

Bedgraph of reads from minus strand (if `make_bedgraph: TRUE`).

Because riboviz aligns to the transcriptome, which represents single-stranded positive-sense RNA, there should be very few reads counted in `minus.bedgraph`. This file is produced by bedtools.


## `plus.bedgraph` 

Bedgraph of reads from plus strand (if `make_bedgraph: TRUE`).

Almost all translated reads should be counted in open reading frames within `plus.bedgraph`, again because riboviz aligns to the transcriptome, which represents single-stranded positive-sense RNA. This file is produced by bedtools.


## `<SAMPLE_ID>.h5` 

Length-sensitive alignments of reads in compressed h5 format. This file is created from the sample Bam file using `Bam_to_H5.R`. Information contained within the h5 file can be accessed using the functions `GetGeneDataMatrix` and `TidyGeneDataMatrix` in R, which will create a tibble showing the number of reads of each length occurring at each position in a gene. More useful functions for working with a h5 file are described in `rscripts/read_count_functions.R`.

## `metagene_start_stop_read_counts.tsv`

A tsv file, showing the sum of reads occurring at positions around the start and stop codons for all genes. Generated by `generate_stats_figs.R` during step “Check for 3nt periodicity globally”, using functions `CalculateThreeNucleotidePeriodicity` and `WriteThreeNucleotidePeriodicity`. These function produce a table showing the number of reads mapping to each position being investigated, using `AllGenes5StartPositionLengthCountsTibble` and `AllGenes3EndPositionLengthCountsTibble` to extract information from the sample h5 file.

```
Pos	Counts	End
-24	536	5'
-23	426	5'
-22	391	5'
-21	428	5'
-20	398	5'
-19	569	5'
```


## `metagene_start_stop_read_counts.pdf`

A meta feature plot showing the total number of reads present within a 75nt window around the start and stop codons. The reads occurring at each position for each gene are summed to give a total for each position, then plotted. It is expected to see a large peak just upstream of the start codon, due to ribosome binding being the slow step of translation. This is typically followed by regular repeating smaller peaks, known as 3nt periodicity, as the majority of reads will map to the first nucleotide of a codon. If the expected features are not seen, then it is possible that there is a problem with annotation files, the adapter listed in the config files or the dataset used is of low quality. 

Good quality 3nt periodicity plot:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/3nt_periodicity.jpg" alt="Good quality 3nt_periodicity plot" width="500"/>

3nt periodicity plot produced using incorrect annotation files:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/bad_3nt_periodicity.jpg" alt="poor quality 3nt_periodicity plot" width="500"/>


## `read_counts_by_length.tsv`

A tsv file showing how many reads are of each length, with the majority being of 28-31nt in length. Generated by `generate_stats_figs.R` during step “Distribution of lengths of all mapped reads” using functions `CalculateReadLengths` and `WriteReadLengths` to get the length and number of reads for each gene, then calculate the total number of reads of each length. These functions extract the information from the h5 file. 

```
Length	Counts
10	0
11	0
12	0
13	0
14	0
15	21567
16	42981
```

## `read_counts_by_length.pdf`

A bar chart showing the lengths of the reads detected in the sample. It is expected that the majority of reads will be 28-31nt long if the adapter sequences have been removed correctly. A good file to check first when running a new dataset, as if the reads peak in the expected range then it is a good indication of a successful run. If no reads are detected then it is a clear indication of something going wrong, such as the wrong UMI expression or adapter sequence being used, leading to reads being unable to align to annotation files.

Example read_counts_by_length plot:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/read_counts_by_length.jpg" alt="poor quality read_counts_by_length plot" width="500"/>


## `pos_sp_nt_freq.tsv`

A tsv file listing the frequency of each different nucleotide at each position of reads of different lengths, including all three possible frames for each read length. Generated by `generate_stats_figs.R` during step “Biases in nucleotide composition along mapped read lengths” using functions `CalculateBiasesInNucleotideComposition` and `WriteBiasesInNucleotideComposition`. For each read length, the frequency of each nucleotide being detected at each position is calculated.

```
Length	Position	Frame	A	C	G	T
10	1	0	0	0	0	0
10	2	0	0	0	0	0
10	3	0	0	0	0	0
10	4	0	0	0	0	0
10	5	0	0	0	0	0
10	6	0	0	0	0	0
```

*Potential names*

read_nt_compostion, nt_freq_read_length, read_nt_freq

## `metagene_normalized_profile_start_stop.pdf`

A plot that shows the mean number of reads mapping to each position upstream and downstream of the start and stop codon for all genes. It is expected to have a peak at the start codon, with the majority of following positions having a relatively consistent mean that is comparatively low to the mean observed at the start codon. 

Example pos_sp_rpf_norm_reads plot:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/pos_sp_rpf_norm_reads.jpg" alt="pos_sp_rpf_norm_reads plot" width="500"/>



## `metagene_normalized_profile_start_stop.tsv`

A tsv file listing the mean and standard deviation of a normalised number of reads detected at positions along genes. The positions are from 1 to 500 nt downstream of the start codon and -500 to -1 nt upstream of the stop codon. Generated by `generate_stats_figs.R` during step “Position specific distribution of reads” using functions `CalculatePositionSpecificDistributionOfReads` and `WritePositionSpecificDistributionOfReads`, which estimate the mean and standard deviation of the number of mapped reads at specific positions of a meta transcript. 

```
Position	Mean	SD	End
1	12.0663784989888	0.395152590437387	5'
2	18.9618100010075	0.574744617743691	5'
3	3.0604586928262		0.130589596966622	5'
4	3.48406026456054	0.131466878328202	5'
5	3.34745502397165	0.160793140986733	5'
6	4.25230808548701	0.146130009657467	5'
```


## `ORF_TPMs_vs_features.tsv`

A tsv file showing the tpm of a set of features and a value for that feature. The features are Length_log10 (log10 of the gene length); uATGs (number of upstream start codons); FE_atg (Free Energy at ATG); FE_cap (Free Energy at the cap); utr (length of the UTR); utr_gc (GC content of the UTR); and polyA (3' polyA tail on mRNA). This is produced by combining the tpms file with the features file, if a features file is provided. This is done by `generate_stats_figs.R` during step "Correlations between TPMs of genes with their sequence-based features" using the functions `CalculateSequenceBasedFeatures` and `WriteSequenceBasedFeatures`.

```
ORF	tpm	Feature	Value
YAL001C	3.42470723290552	Length_log10	3.06445798922692
YAL002W	1.78650805583587	Length_log10	3.10516942799933
YAL003W	4043.513688639	  Length_log10	2.31386722036915
YAL007C	105.576584188032	Length_log10	2.33243845991561
YAL008W	11.8454807087593	Length_log10	2.29666519026153
YAL010C	4.82007633057095	Length_log10	2.69284691927723
```

Only output if `--features-file` was defined.


## `ORF_TPMs_vs_features.pdf` 

The features pdf relates the tpm value of different genes to a variety of different sequence features. This highlights any trends in feature value as the tpm changes. Values for the different features come from the features-file, if it is provided and are described above.

Example features plot:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/features.jpg" alt=" features plot" width="500"/>

Only output if `--features-file` was defined.

 

## `tpms.tsv`

A tsv file listing the rpb (reads per base) and tpm of the ORFs of a sample, along with the number of reads detected in the ORF. Generated by `generate_stats_figs.R` during the step “Calculate TPMs of genes” using functions `CalculateGeneTranscriptsPerMillion` and `WriteGeneTranscriptsPerMillion`. These functions get the total number of reads per gene from the h5 file which they use to calculate the rpb and tpm.

```
ORF	readcount	rpb	tpm
Q0045	17	0.0102905569007264	0.665266152991222
Q0050	28	0.0109717868338558	0.709306453364419
Q0055	6	0.00229709035222052	0.148502794983537
Q0060	3	0.00231660231660232	0.1497642086861
Q0065	26	0.0151338766006985	0.97837813474214
Q0070	20	0.0103092783505155	0.66647646133505
```

*Potential names*

<SAMPLE_ID>_tmps

## `normalized_density_APEsites_per_codon.tsv` 

A tsv file showing the correlatation of different codons to features based on a provided tRNA file, which gives the Amino acids, the tRNA estimates, the tAI (tRNA Adaptation Index), Microarray values, and RNA.seq values for each codon. These are used to calculate mean ribosome-densities at the A/P/E sites for each codon.

Produced by `generate_stats_figs` during step "Codon-specific ribosome densities for correlations with tRNAs" using functions `CalculateCodonSpecificRibosomeDensityTRNACorrelation` and `WriteCodonSpecificRibosomeDensityTRNACorrelation`.

```
AA	Codon	tRNA	tAI	Microarray	RNA.seq	A	P	E
K	AAA	7	0.431034	222273	82386	0.701009246604793	0.757271500717756	0.882163216975167
N	AAC	10	0.615764	378101	110849	0.608843172114343	0.868064826175943	0.844584176727076
K	AAG	14	1	        397111	83036	0.719699425589095	0.985491091356775	1.11948381450794
N	AAT	6.4	0.27032	        241984.64 70943.36	0.566975111040301	0.786255203154437	0.76039088714112
T	ACA	4	0.246373	105862	47598	0.933640392583143	1.007034200634	1.11922337626619
```

Only output if `--t-rna-file` and `--codon-positions-file` were defined.


## `normalized_density_APEsites_per_codon_long.tsv`

A tsv file showing the correlatation of different codons to features based on a provided tRNA file, which gives the Amino acids, the tRNA type (tRNA, tAI, Microarray, and RNA.seq), the tRNA values, the site in the ribosome and the ribodensity for each codon. It is produced by reformatting the normalized_density_APEsites_per_codon.tsv file. 

```
AA	Codon	tRNA_type	tRNA_value	Site	Ribodens
K	AAA	tRNA	7	A	1.19938409239398
K	AAA	tRNA	7	P	1.2109078995196
K	AAA	tRNA	7	E	1.02770556029175
K	AAA	tAI	0.431034	A	1.19938409239398
K	AAA	tAI	0.431034	P	1.2109078995196
K	AAA	tAI	0.431034	E	1.02770556029175
K	AAA	Microarray	222273	A	1.19938409239398
K	AAA	Microarray	222273	P	1.2109078995196
K	AAA	Microarray	222273	E	1.02770556029175
K	AAA	RNA.seq	82386	A	1.19938409239398
K	AAA	RNA.seq	82386	P	1.2109078995196
K	AAA	RNA.seq	82386	E	1.02770556029175

```

Only output if `--t-rna-file` and `--codon-positions-file` were defined.



## `normalized_density_APEsites_per_codon.pdf` 

This plot shows a range of features and relates them to ribosome densitiy on the A, P and E sites. 4 features are shown; Microarry, RNA.seq, tAI and tRNA. These are taken from the `--t-rna-file` if provided. Each codon codon has a different value for each of these features, described generally here as a tRNA_value. These tRNA_values are plotted against the ribodensity at each site of the ribosome, showing any relationships or trends.

Example normalized_density_apesites_per_codon plot:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/codon_ribodens.jpg" alt="codon_ribodens plot" width="500"/>

Only output if `--t-rna-file` and `--codon-positions-file` were defined.


## `gene_position_length_counts_5start.tsv`

A TSV file containing the number of reads of different lengths mapping to positions around the start codon. This file is created by `generate_stats_figs.R` during step "Check for 3nt periodicity globally" by functions CalculateGenePositionLengthCounts5Start and WriteGenePositionLengthCounts5Start. 

```
ReadLen	Pos	Counts
10	-24	0
11	-24	0
12	-24	0
13	-24	0
14	-24	0
15	-24	3
16	-24	4
17	-24	9

```

## `metagene_start_barplot_by_length.pdf`

A meta-feature bar chart showing the number of reads occurring at positions around the start codons of genes, faceted by read length. As the majority of reads will be 28-31nt in length, only bar charts for lengths 26-32nt are shown. It is expected that each length will show a peak of reads just upstream of the start codon for all lengths, and then an observable 3nt periodicity following the peak, which will be more distinct in the more common read lengths. This figure is created as part of `generate_stats_figs.R` step “Check for 3nt periodicity globally”, using data that is saved in the `3nt_periodicity.tsv` file. Note: the Y axis scale will vary for each read length. 

Example metagene_start_barplot_by_length.pdf:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/startcodon_ribogridbar.jpg" alt="startcodon_ribogridbar plot" width="500"/>


## `metagene_start_ribogrid_by_length.pdf`

A meta-feature heatmap showing the number of reads occurring at positions around the start codons of genes with the y axis as read length, and the colour intensity showing read count. It is expected that each length will show a peak of reads just upstream of the start codon for all lengths, shown by an intense dark purple, and then an observable 3nt periodicity following the peak, which will be more distinct in the more common read lengths. This figure is created as part of `generate_stats_figs.R` step “Check for 3nt periodicity globally”, using data that is saved in the `3nt_periodicity.tsv` file. 

Example metagene_start_ribogrid_by_length.pdf:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/startcodon_ribogrid.jpg" alt="startcodon_ribogrid plot" width="500"/>


## `read_frame_per_ORF.tsv` 

A tsv file showing the count threshold for each frame and the p value comparing the counts in frame 0 to frames 1, 2 and both 1 and 2 for each gene. It is expected that there will be a significant difference between the frames, with the majority of reads mapping to frame 0. Generated by `generate_stats_figs.R` during step “Check for 3nt periodicity (frame) by Gene” using functions `CalculateGeneReadFrames` and `WriteGeneReadFrames`. These functions count the number of reads in each frame and calculates the Wilcoxon rank-sum paired test p-values.

```
gene	Ct_fr0	Ct_fr1	Ct_fr2	pval_fr0vs1	pval_fr0vs2	pval_fr0vsboth
YAL068C	0	0	0	1	1	1
YAL067W-A	0	0	0	1	1	1
YAL067C	1	0	0	0.5	0.5	0.5
YAL065C	65	10	56	0.000225406271329839	0.433008955500838	0.00465962908671896
YAL064W-B	0	0	0	1	1	1
YAL064C-A	0	0	0	1	1	1
```

Only output if `--asite-disp-length-file` was defined.

## `read_frame_per_ORF_filtered.tsv`

A tsv file showing the count threshold for each frame and the p value comparing the counts in frame 0 to frames 1, 2 and both 1 and 2 for each gene. It is expected that there will be a significant difference between the frames, with the majority of reads mapping to frame 0. This file is produced by filtering the read_frame_per_ORF.tsv file, keeping only the genes where Ct_fr0 + Ct_fr1 + Ct_fr2 is greater than the count_threshold, which is defined in the config file.

```
gene	Ct_fr0	Ct_fr1	Ct_fr2	pval_fr0vs1	pval_fr0vs2	pval_fr0vsboth
YAL065C	65	10	56	0.000225406271329839	0.433008955500838	0.00465962908671896
YAL063C	859	207	944	2.71355108503869e-21	0.976833449475798	8.09866408105508e-08
YAL060W	267	40	201	7.93480334492926e-15	0.0274846446691625	5.69179781323741e-12
YAL059W	173	34	159	2.29896520195999e-09	0.339582970614079	4.5634250757101e-06
YAL058W	43	4	31	3.99564191533441e-08	0.175807992651532	4.21595642555769e-06
```

Only output if `--asite-disp-length-file` was defined.

## `frame_proportions_per_ORF.pdf` 

A box plot of the proportion of reads mapping to each nucleotide of a codon for each ORF, with the first codon being denoted as Frame 0, the second being Frame 1 and the 3rd being Frame 2. It is expected that the majority of reads will be mapping to the first nucleotide, Frame 0, so this box will be higher than the others. 

Example 3ntframe_propbygene.pdf:

<img src="https://github.com/3mma-mack/Riboviz-honours/blob/main/riboviz_images/3ntframe_propbygene.jpg" alt="3ntframe_propbygene plot" width="500"/>

Only output if `--asite-disp-length-file` was defined.
