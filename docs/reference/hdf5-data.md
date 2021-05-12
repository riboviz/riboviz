# Structure of HDF5 data

For analyses of ribosome footprinting and RNA-seq datasets, we store summaries of aligned read data in Hierarchical Data Format ([HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)) format. HDF5 allows for rapid access to mapped reads of a particular length to any coding sequence.

To learn more about accessing and manipulating HDF5 files in R, see Bioconductor's [HDF5 interface to R](https://bioconductor.org/packages/release/bioc/html/rhdf5.html) and Neon's [Introduction to HDF5 Files in R](https://www.neonscience.org/hdf5-intro-r).

---

## `bam_to_h5.R`

Given a GFF file and a BAM file, `bam_to_h5.R` creates an HDF5 file with information about a feature (e.g. CDS, ORF, or uORF).

It is passed the following configuration parameters from the RiboViz configuration (see [Configuring the RiboViz workflow](../user/prep-riboviz-config.md)):

| Parameter | Description |
| --------- | ----------- |
| `buffer` | Length of flanking region around the feature |
| `dataset` | Human-readable name of the dataset |
| `is_riboviz_gff` | Does `orf_gff_file` contain 3 elements per gene - UTR5, feature, and UTR3? |
| `max_read_length` | Maximum read length in H5 output |
| `min_read_length` | Minimum read length in H5 output |
| `num_processes` | Number of processes to parallelize over, used by specific steps in the workflow |
| `orf_gff_file` | Matched genome feature file, specifying coding sequences locations (start and stop coordinates) within the transcripts (GTF/GFF3 file) |
| `primary_id` | Primary gene IDs to access the data (YAL001C, YAL003W, etc.) |
| `secondary_id` | Secondary gene IDs to access the data (COX1, EFB1, etc. or `NULL`) |
| `stop_in_cds` | Are stop codons part of the feature annotations in `orf_gff_file`? If both `stop_in_feature` and `stop_in_cds` are defined then `stop_in_feature` takes precedence. |
| `stop_in_feature` | Are stop codons part of the feature annotations in `orf_gff_file`? If not provided and `stop_in_cds` is provided then the value of `stop_in_cds` is used for `stop_in_feature`. If both `stop_in_feature` and `stop_in_cds` are defined then `stop_in_feature` takes precedence. |

At present, the default feature is assumed to be `CDS`.

All reads are mapped to their 5' ends.

`primary_id` is the name of an attribute in `orf_gff_file` (e.g. `Name`) expected to hold gene names.

`secondary_id` is the name of an attribute in `orf_gff_file` (e.g. `ID`) expected to hold gene names or, if none, `NULL`. If provided then these alternative gene names are used to create symbolic links in the H5 files to the entries for each  gene.

`primary_id` is used as GFF column e.g. "Name".

`secondary_id` is used as a GFF column to get alternative gene names which are associated with gene names. These are used to create symbolic links in the H5 file to the entries for the original genes.

If `is_riboviz_gff` then:

* Feature (e.g. `CDS`, `ORF`, or `uORF`), `UTR5` and `UTR3` entries from `orf_gff_file` are used.
* `buffer` is ignored.

If not `is_riboviz_gff` then:

* Feature (e.g. `CDS`, `ORF`, or `uORF`) entries from `orf_gff_file` are used.
* `UTR5` and `UTR3` entries from `orf_gff_file` are ignored.
* `buffer` is used as the width of left and right flanks.
* `stop_in_feature` states where the stop codon is located (true if within the feature, false otherwise).

---

## HDF5 file structure

The `reads` group, `/<gene>/<dataset>/reads`, for a `<gene>` has several attributes associated with it. These are summary statistics and other information about the gene and dataset within the `reads` group. The list of attributes are as follows.

| Attribute | Description | Origin |
| --------- |------------ | ------ |
| `buffer_left` | Number of nucleotides upstream of the start codon (ATG) (UTR5 length) | EITHER position of start codon (from `orf_gff_file`) - 1 OR, if `is_riboviz_gff` is false, `buffer` |
| `buffer_right` | Number of nucleotides downstream of the stop codon (TAA/TAG/TGA) (UTR3 length) | From `orf_gff_file` OR, if `is_riboviz_gff` is false, `buffer` + feature length from `orf_gff_file` + `buffer` - `stop_codon_pos[3]` (see below) |
| `start_codon_pos` | Positions corresponding to start codon of feature (typically 251, 252, 253) | From `orf_gff_file` |
| `stop_codon_pos` | Positions corresponding to stop codon of feature | Positions from `orf_gff_file`. If `is_riboviz_gff` is false and `stop_in_feature` is false, position is last nucleotide of feature + 1 |
| `lengths` | Lengths of mapped reads | Range from `min_read_length` to `max_read_length` (e.g. 10,...,50) |
| `reads_by_len` | Counts of number of ribosome sequences of each length | BAM file |
| `reads_total` | Total number of ribosome sequences | BAM file. Equal to number of non-zero reads in `reads_by_len` |

The HDF5 file is organized in a hierarchy, `/<gene>/<dataset>/reads/data`, where `<gene>` is each gene name in the BAM file and `orf_gff_file` passed to `bam_to_h5.R`.

The `data` table in `/<gene>/<dataset>/reads/data`, for a `<gene>` has the positions and lengths of ribosome sequences within the organism data (determined from the BAM file). It is an integer table with each row representing a read length and columns representing nucleotide positions. The first row corresponds to reads of length `min_read_length` and the last row corresponds to reads of length `max_read_length`.

A template HDF5 file, showing how the HDF5 file relates to information in the RiboViz configuration (and `bam_to_h5.R command-line parameters), `orf_gff_fasta` and a BAM file is as follows:

```
HDF5 "<file>.h5" {
GROUP "/" {
   GROUP "<gene>" {
      GROUP "<dataset>" {
	 GROUP "reads" {
	    ATTRIBUTE "buffer_left" {
	       DATATYPE  H5T_STD_I32LE
	       DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
	       DATA {
	       (0): <Feature UTR5 length from orf_gff_file OR buffer (if is_riboviz_gff is false)>
	       }
	    }
	    ATTRIBUTE "buffer_right" {
	       DATATYPE  H5T_STD_I32LE
	       DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
	       DATA {
	       (0): <Feature UTR3 length from orf_gff_file OR buffer (if is_riboviz_gff is false)>
	       }
	    }
	    ATTRIBUTE "lengths" {
               DATATYPE  H5T_STD_I32LE
               DATASPACE  SIMPLE { (<read_length>) / (<read_length>) }
               DATA {
               (0): <min_read_length>, <min_read_length> + 1, ... , <min_read_length> + <m> - 1,
               (<m>): <min_read_length>+<m>, <min_read_length> + <m> + 1, ... , <min_read_length> + <n> - 1,
               (<n>): <min_read_length>+<n>, <min_read_length> + <n> + 1, ... , <min_read_length> + <max_read_length>
               }
            }
            ATTRIBUTE "reads_by_len" {
               DATATYPE  H5T_STD_I32LE
               DATASPACE  SIMPLE { (<read_length>) / (<read_length>) }
               DATA {
               (0): <see below>
               (<m>): <see below>
               (<n>): <see below>
               }
            }
           ATTRIBUTE "reads_total" {
               DATATYPE  H5T_STD_I32LE
               DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
               DATA {
               (0): <number of non-zero values in reads_by_len>
               }
            }
            ATTRIBUTE "start_codon_pos" {
               DATATYPE  H5T_STD_I32LE
               DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
               DATA {
               (0): <position of 1st nt of feature start codon from orf_gff_file>,
		    <position of 2nd nt of feature start codon from orf_gff_file>,
		    <position of 3rd nt of feature start codon from orf_gff_file>
               }
            }
            ATTRIBUTE "stop_codon_pos" {
               DATATYPE  H5T_STD_I32LE
               DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
               DATA {
               (0): <position of 1st nt of feature stop codon from orf_gff_file>,
	            <position of 2nd nt of feature stop codon from orf_gff_file>,
	            <position of 3rd nt of feature stop codon from orf_gff_file>
               }
            }
            DATASET "data" {
               DATATYPE  H5T_STD_I32LE
               DATASPACE  SIMPLE { ( <sequence_length>, <read_length> ) /
	                           ( <sequence_length>, <read_length> ) }
               DATA {
               (0, 0): <see below>
               (0, <m>): <see below>
               (0, <n>): <see below>
               ...
               (<p>, 0): <see below>
               (<p>, <m>): <see below>
               (<p>, <n>): <see below>
               ...
               (<sequence_length - 1>, 0): <see below>
               (<sequence_length - 1>, <m>): <see below>
               (<sequence_length - 1>, <n>): <see below>
               }
            }
         }
      }
   }
   ...
}
```

where:


* `gene`: Gene ID from `orf_gff_file` `primary_id` attribute.
* `read_length`: `max_read_length` - `min_read_length` + 1
* 0 < `m` < `n` < `read_length`. Exact values of `m` and `n` may differ across specific `ATTRIBUTE` and `DATA` items.
* `reads_by_len`:
  - `reads_by_len[i]` = number of alignments in BAM file which have `Flag` equal to 0 or 256 and length equal to `lengths[i]`. This equals the sum of `DATA[*, i]` i.e. sum across all positions for a specific read length.
* `sequence_length`: position of final nucleotide of UTR3 from `orf_gff_file` (equal to length of sequence from BAM file header `LN` value).
  - If `is_riboviz_gff` is false when this is equal to `buffer` + feature length + `buffer`.
* `DATASET "data"`:
  - 0 <= `p` <= `sequence_length - 1`
  - `DATA[p, i]` equals 1 if there is an alignment in the BAM file at position `p`+1 which has length equal to `lengths[i]` and alignment has `Flag` value 0 or 256 (see [Understanding the BAM flags](https://davetang.org/muse/2014/03/06/understanding-bam-flags/) and [Decoding SAM flags](https://broadinstitute.github.io/picard/explain-flags.html)); 0 otherwise.

---

## Example

A snippet from an example HDF5 file is shown below. 

```
                              group              name       otype  dclass       dim
0                                 /           YAL001C   H5I_GROUP                  
1                          /YAL001C 2016_Weinberg_RPF   H5I_GROUP                  
2        /YAL001C/2016_Weinberg_RPF             reads   H5I_GROUP                  
3  /YAL001C/2016_Weinberg_RPF/reads              data H5I_DATASET INTEGER 36 x 3980
4                                 /           YAL002W   H5I_GROUP                  
5                          /YAL002W 2016_Weinberg_RPF   H5I_GROUP                  
6        /YAL002W/2016_Weinberg_RPF             reads   H5I_GROUP                  
7  /YAL002W/2016_Weinberg_RPF/reads              data H5I_DATASET INTEGER 36 x 4322
8                                 /           YAL003W   H5I_GROUP                  
9                          /YAL003W 2016_Weinberg_RPF   H5I_GROUP                  
10       /YAL003W/2016_Weinberg_RPF             reads   H5I_GROUP                  
11 /YAL003W/2016_Weinberg_RPF/reads              data H5I_DATASET INTEGER 36 x 1118
```

The figure below shows an example where `min_read_length` is 15 and `max_read_length` is 50:

![H5 architecture](../images/h5_architecture.png)
