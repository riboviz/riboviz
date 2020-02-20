# Structure of HDF5 data

For all subsequent analyses of ribosome footprinting and RNA-seq datasets, we store summaries of aligned read data in Hierarchical Data Format ([HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)) format. HDF5 allows for rapid access to mapped reads of a particular length to any coding sequence.

To learn more about accessing and manipulating HDF5 files in R, see Bioconductor's [HDF5 interface to R](https://bioconductor.org/packages/release/bioc/html/rhdf5.html) and Neon's [Introduction to HDF5 Files in R](https://www.neonscience.org/hdf5-intro-r).

## HDF5 file architecture

The HDF5 file is organized in the following hierarchy `/<Gene>/<Dataset>/reads/data`. A snippet from an example HDF5 file is shown below. 

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
The `data` table is an integer table with each rows representing a read length and columns representing nucleotide positions. The first row corresponds to reads of length 15 and the last row corresponds to reads of length 50. All reads are mapped to their 5' ends (see below).

![as](https://gdurl.com/9UUu)

## Read attributes

The `reads` group in `/<Gene>/<Dataset>/reads/data` have several attributes associated with it. These are summary statistics and other information about the gene and dataset within the `reads` group. The list of attributes are

* `reads_total` : Sum of reads of all lenghts between -25 to +25 of a CDS
* `buffer_left` : Number of nucleotides upstream of the start codon (ATG) - 250nt
* `buffer_right` : Number of nucleotides downstream of the stop codon (TAA/TAG/TGA) - 247nt
* `start_codon_pos` : Positions corresponding to the start codon - (251,252,253)
* `stop_codon_pos` : Positions corresponding to the stop codon (variable)
* `reads_by_len` : Sum of reads between -25 to +25 of a CDS for each length
* `lengths` : Lengths of mapped reads (15-50)
