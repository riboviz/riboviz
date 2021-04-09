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
| `is_riboviz_gff` | Does the GFF file contain 3 elements per gene - UTR5, feature, and UTR3? |
| `max_read_length` | Maximum read length in H5 output |
| `min_read_length` | Minimum read length in H5 output |
| `num_processes` | Number of processes to parallelize over, used by specific steps in the workflow |
| `orf_gff_file` | Matched genome feature file, specifying coding sequences locations (start and stop coordinates) within the transcripts (GTF/GFF3 file) |
| `primary_id` | Primary gene IDs to access the data (YAL001C, YAL003W, etc.) |
| `secondary_id` | Secondary gene IDs to access the data (COX1, EFB1, etc. or `NULL`) |
| `stop_in_cds` | Are stop codons part of the feature annotations in GFF? |

At present, the default feature is assumed to be `CDS`.

All reads are mapped to their 5' ends.

`primary_id` is used as GFF column e.g. "Name".

`secondary_id` is used as a GFF column to get alternative gene names which are associated with gene names. These are used to create symbolic links in the H5 file to the entries for the original genes.

If `is_riboviz_gff == TRUE` then:

* Feature (e.g. `CDS`, `ORF`, or `uORF`), `UTR5` and `UTR3` entries from the GFF are used.
* `buffer` is ignored.

If `is_riboviz_gff == FALSE` then:

* Feature (e.g. `CDS`, `ORF`, or `uORF`), `UTR5` and `UTR3` entries from the GFF are used.
* `UTR5` and `UTR3` entries from the GFF are ignored.
* `buffer` is used as the width of left and right flanks, and to calculate start and stop codon locations for genes.
* Whether `stop_in_cds` is `TRUE` or `FALSE` also affects the calculations for stop codon locations.

---

## HDF5 file structure

The HDF5 file is organized in the following hierarchy `/<gene>/<dataset>/reads/data`, where `<gene>` is each gene name in the GFF and BAM files passed to `bam_to_h5.R`.

The `reads` group, `/<gene>/<dataset>/reads`, for a `<gene>` has several attributes associated with it. These are summary statistics and other information about the gene and dataset within the `reads` group. The list of attributes are as follows.

| Attribute | Description | Origin |
| --------- |------------ | ------ |
| `buffer_left` | Number of nucleotides upstream of the start codon (ATG) (UTR5 length) | EITHER position of start codon (from GFF file) - 1 OR, if `is_riboviz_gff` is false, `buffer` |
| `buffer_right` | Number of nucleotides downstream of the stop codon (TAA/TAG/TGA) (UTR3 length) | GFF file OR, if `is_riboviz_gff` is false, `buffer` + feature length from GFF file + `buffer` - `stop_codon_pos[3]` (see below) |
| `start_codon_pos` | Positions corresponding to start codon of feature (typically 251, 252, 253) | EITHER positions from GFF file OR, if `is_riboviz_gff` is false, `buffer+1,...,buffer+3` |
| `stop_codon_pos` | Positions corresponding to stop codon of feature | EITHER positions from GFF file OR, if `is_riboviz_gff` is false, `p,...,p+2` where `p` is `buffer` + feature length from GFF file - `offset` where `offset` is 2 if `stop_in_cds` is true or -1 otherwise |
| `lengths` | Lengths of mapped reads | Range from `min_read_length` to `max_read_length` (e.g. 10,...,50) |
| `reads_by_len` | Counts of number of ribosome sequences of each length | BAM file |
| `reads_total` | Total number of ribosome sequences | BAM file. Equal to number of non-zero reads in `reads_by_len` |

The `data` table in `/<gene>/<dataset>/reads/data`, for a `<gene>` has the positions and lengths of ribosome sequences within the organism data (determined from the BAM file). It is an integer table with each row representing a read length and columns representing nucleotide positions. The first row corresponds to reads of length `min_read_length` and the last row corresponds to reads of length `max_read_length`.

A template HDF5 file, showing how the HDF5 file relates to information in the RiboViz configuration (and `bam_to_h5.R command-line parameters), GFF and BAM files is as follows:

```
HDF5 "<file>.h5" {
GROUP "/" {
   GROUP "<gene>" {
      GROUP "<dataset>" {
	 GROUP "reads" {
	    ATTRIBUTE "buffer_left" {
	       DATATYPE  H5T_IEEE_F64LE
	       DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
	       DATA {
	       (0): <Feature UTR5 length from GFF file OR buffer (if is_riboviz_gff is false)>
	       }
	    }
	    ATTRIBUTE "buffer_right" {
	       DATATYPE  H5T_STD_I32LE
	       DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
	       DATA {
	       (0): <Feature UTR3 length from GFF file OR buffer (if is_riboviz_gff is false)>
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
               DATATYPE  H5T_IEEE_F64LE
               DATASPACE  SIMPLE { (<read_length>) / (<read_length>) }
               DATA {
               (0): <see below>
               (<m>): <see below>
               (<n>): <see below>
               }
            }
           ATTRIBUTE "reads_total" {
               DATATYPE  H5T_IEEE_F64LE
               DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
               DATA {
               (0): <number of non-zero values in reads_by_len>
               }
            }
            ATTRIBUTE "start_codon_pos" {
               DATATYPE  H5T_STD_I32LE
               DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
               DATA {
               (0): <position of 1st nt of feature start codon from GFF file>,
		    <position of 2nd nt of feature start codon from GFF file>,
		    <position of 3rd nt of feature start codon from GFF file>
               }
            }
            ATTRIBUTE "stop_codon_pos" {
               DATATYPE  H5T_STD_I32LE
               DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
               DATA {
               (0): <position of 1st nt of feature stop codon from GFF file>,
	            <position of 2nd nt of feature stop codon from GFF file>,
	            <position of 3rd nt of feature stop codon from GFF file>
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

* `gene`: Gene ID from GFF file (equal to reference sequence name from BAM file).
* `read_length`: `max_read_length` - `min_read_length` + 1
* 0 < `m` < `n` < `read_length`. Exact values of `m` and `n` may differ across specific `ATTRIBUTE` and `DATA` items.
* `reads_by_len`:
  - `reads_by_len[i]` = number of alignments in BAM which have Flag equal to 0 or 256 and length equal to `lengths[i]`. This equals the sum of `DATA[*, i]` i.e. sum across all positions for a specific read length.
* `sequence_length`: position of final nt of UTR3 from GFF file (equal to length of sequence from BAM file header `LN` value).
  - If `is_riboviz_gff` is `false` when this is equal to `buffer` + feature length + `buffer`.
* `DATASET "data"`:
  - 0 <= `p` <= `sequence_length - 1`
  - `DATA[p, i]` equals 1 if there is an alignment in the BAM file at position `p`+1 which has length equal to `lengths[i]` and BAM alignment has Flag 0 or 256 (see [Understanding the BAM flags](https://davetang.org/muse/2014/03/06/understanding-bam-flags/) and [Decoding SAM flags](https://broadinstitute.github.io/picard/explain-flags.html)); 0 otherwise.

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
