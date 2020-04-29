# load library {here}
suppressMessages(library(here, quietly=T))

# source generate_stats_figs.R
source(here::here("rscripts", "generate_stats_figs.R"))

print("testing h5-open-106 fixes")

#### GetGeneDataMatrix() ####

print("testing GetGeneDataMatrix()")
# # function to get data matrix of read counts for gene and dataset from hdf5file


# OLD VERSION:
# hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
#
# GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
#   hdf5file %>%
#     rhdf5::H5Dopen(
#       name = paste0("/", gene, "/", dataset, "/reads/data")
#     ) %>%
#     rhdf5::H5Dread() %>%
#     return()
# }


# NEW VERSION
# hdf5file <- hd_file # filehandle for the h5 file
# 
# GetGeneDatamatrixAlt <- function(gene, dataset, hdf5file){
#   rhdf5::h5read(file = hdf5file, name = paste0("/", gene, "/", dataset, "/reads/data")) %>%
#     return()
# }


# old version outputs:
GeneDatamatrixYAL001C_OLD <- GetGeneDatamatrix("YAL001C", "vignette", hdf5file)

str(GeneDatamatrixYAL001C_OLD)
# int [1:41, 1:3983] 0 0 0 0 0 0 0 0 0 0 ...

# close connections to avoid confusion
h5closeAll()


# RUN NEW VERSION CODE:
hdf5file <- hd_file # filehandle for the h5 file

GetGeneDatamatrixAlt <- function(gene, dataset, hdf5file){
  rhdf5::h5read(file = hdf5file, name = paste0("/", gene, "/", dataset, "/reads/data")) %>%
    return()
}

# run new version
GeneDatamatrixYAL001C_NEW <- GetGeneDatamatrixAlt("YAL001C", "vignette", hdf5file)

# check output structure
str(GeneDatamatrixYAL001C_NEW)
# int [1:41, 1:3983] 0 0 0 0 0 0 0 0 0 0 ...

# close connections to avoid confusion
h5closeAll()


# are outputs equal between old and new?
all.equal(GeneDatamatrixYAL001C_OLD, GeneDatamatrixYAL001C_NEW)
# [1] TRUE

print("GetGeneDataMatrix() old and new equal? ")
print(all.equal(GeneDatamatrixYAL001C_OLD, GeneDatamatrixYAL001C_NEW))



#### gene_sp_read_length ####

print("testing gene_sp_read_length")

# So this section equates to the function GetGeneReadLength() in 
# read_counts_functions.R in branch regen-stats-32, and is because I haven't 
# changed out the rhdf5::H5Aread(), rhdf5::H5Aopen(), rhdf5::H5Gopen()` functions yet.
# in branch regen-stats-32 this is currently:
# GetGeneReadLength <- function(gene, dataset, hdf5file) {
#   hdf5file %>%
#     rhdf5::H5Gopen(
#       name = paste0("/", gene, "/", dataset, "/reads")
#     ) %>%
#     rhdf5::H5Aopen(name = "reads_by_len") %>%
#     rhdf5::H5Aread() %>% # return array with the read data
#     return()
# }


# OLD VERSION CODE:
# hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
#
# gene_sp_read_length <- lapply(gene_names, function(x) {
#   +   rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", x, "/", dataset, "/reads")), "reads_by_len"))
#   + })


# NEW VERSION CODE: 
# hdf5file <- hd_file # filehandle for the h5 file
# 
# GetGeneReadLength <- function(gene, hdf5file){ # new version of regen_stats_32 branch function in read_count_functions.R
#   rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["reads_by_len"]]
# }
# 
# gene_sp_read_length_Alt <- lapply(gene_names, function(x) {
#   GetGeneReadLength(x, hdf5file = hdf5file)
# })
  

# RUN OLD VERSION: 
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
gene_sp_read_length <- lapply(gene_names, function(x) {
   rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", x, "/", dataset, "/reads")), "reads_by_len"))
})

# close connections to avoid confusion
h5closeAll()

# check output structure etc
# str(gene_sp_read_length)
# List of 68
# $ : num [1:41(1d)] 0 0 0 0 0 1 0 0 0 0 ...
# $ : num [1:41(1d)] 0 0 0 0 0 0 0 0 0 0 ...
# $ : num [1:41(1d)] 0 0 0 0 1 0 3 3 9 8 ...
# $ : num [1:41(1d)] 0 0 0 3 5 3 15 16 17 29 ...
# ...
length(gene_sp_read_length)
# [1] 68


# RUN NEW VERSION
hdf5file <- hd_file # filehandle for the h5 file

GetGeneReadLength <- function(gene, hdf5file){
  rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["reads_by_len"]]
}

gene_sp_read_length_Alt <- lapply(gene_names, function(x) {
  GetGeneReadLength(x, hdf5file = hdf5file)
})

# check output structure etc
# str(gene_sp_read_length_Alt)
# List of 68
# $ : num [1:41(1d)] 0 0 0 0 0 1 0 0 0 0 ...
# $ : num [1:41(1d)] 0 0 0 0 0 0 0 0 0 0 ...
# $ : num [1:41(1d)] 0 0 0 0 1 0 3 3 9 8 ...
# $ : num [1:41(1d)] 0 0 0 3 5 3 15 16 17 29 ...
# ...
length(gene_sp_read_length_Alt)
# [1] 68

# close connections to avoid confusion
h5closeAll()


# COMPARE OLD AND NEW:
all.equal(gene_sp_read_length, gene_sp_read_length_Alt)

print("gene_sp_read_length old and new equal? ")
print(all.equal(gene_sp_read_length, gene_sp_read_length_Alt))



#### GetGeneLength() ####

print("testing GetGeneLength()")


# OLD VERSION CODE: 
# hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# # read length-specific total counts stored as attributes of 'reads_total' in H5 file
# GetGeneLength <- function(gene, hdf5file) {
#   start_codon_pos <- rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "start_codon_pos"))[1]
#   stop_codon_pos <- rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "stop_codon_pos"))[1]
#   return(stop_codon_pos - start_codon_pos)
# }


# NEW VERSION CODE: 
# hdf5file <- hd_file # filehandle for the h5 file
# GetGeneLengthAlt <- function(gene, hdf5file) {
#   start_codon_pos <- rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["start_codon_pos"]][1]
#   stop_codon_pos <- rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["stop_codon_pos"]][1]
#   return(stop_codon_pos - start_codon_pos)
# }


# RUN OLD VERSION:
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
GetGeneLength <- function(gene, hdf5file) {
  start_codon_pos <- rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "start_codon_pos"))[1]
  stop_codon_pos <- rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "stop_codon_pos"))[1]
  return(stop_codon_pos - start_codon_pos)
}
GeneLengthYAL001C_OLD <- GetGeneLength("YAL001C", hdf5file)

# close connections to avoid confusion
h5closeAll()

# check output structure etc
str(GeneLengthYAL001C_OLD)
# int 3480


# RUN NEW VERSION CODE: 
hdf5file <- hd_file # filehandle for the h5 file
GetGeneLengthAlt <- function(gene, hdf5file) {
  start_codon_pos <- rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["start_codon_pos"]][1]
  stop_codon_pos <- rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["stop_codon_pos"]][1]
  return(stop_codon_pos - start_codon_pos)
}
GeneLengthYAL001C_NEW <- GetGeneLengthAlt("YAL001C", hdf5file)

# close connections to avoid confusion
h5closeAll()

# check output structure etc
str(GeneLengthYAL001C_NEW)
# int 3480


# COMPARE OLD AND NEW:
all.equal(GeneLengthYAL001C_OLD, GeneLengthYAL001C_NEW)

print("GetGeneLength() old and new equal? ")
print(all.equal(GeneLengthYAL001C_OLD, GeneLengthYAL001C_NEW))



#### GetGeneReadsTotal() ####

print("testing GetGeneReadsTotal()")


# OLD VERSION CODE:
# hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# GetGeneReadsTotal <- function(gene, hdf5file) {
#   rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "reads_total"))
# }


# NEW VERSION CODE:
# hdf5file <- hd_file # filehandle for the h5 file
# GetGeneReadsTotalAlt <- function(gene, hdf5file){
#   rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[5]
# }


# RUN OLD VERSION:
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
GetGeneReadsTotal <- function(gene, hdf5file) {
  rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "reads_total"))
}
GeneReadsTotalYAL001C_OLD <- GetGeneReadsTotal("YAL001C", hdf5file)

# close connections to avoid confusion
h5closeAll()

# check output structure etc
str(GeneReadsTotalYAL001C_OLD)
# num [1(1d)] 4
head(GeneReadsTotalYAL001C_OLD)
# [1] 4
length(GeneReadsTotalYAL001C_OLD)
#[1] 1


# RUN NEW VERSION CODE:
hdf5file <- hd_file # filehandle for the h5 file
GetGeneReadsTotalAlt <- function(gene, hdf5file){
  rhdf5::h5readAttributes(hd_file, name=paste0("/", gene, "/", dataset, "/reads"))[["reads_total"]]
}
GeneReadsTotalYAL001C_NEW <- GetGeneReadsTotalAlt("YAL001C", hdf5file)

# close connections to avoid confusion
h5closeAll()

# check output structure etc
str(GeneReadsTotalYAL001C_NEW)
# num [1(1d)] 4

head(GeneReadsTotalYAL001C_NEW)
# [1] 4

length(GeneReadsTotalYAL001C_NEW)
# [1] 1


# are outputs equal between old and new?
all.equal(GeneReadsTotalYAL001C_OLD, GeneReadsTotalYAL001C_NEW)

print("GetGeneReadsTotal() old and new equal? ")
print(all.equal(GeneReadsTotalYAL001C_OLD, GeneReadsTotalYAL001C_NEW))



###

print("completed testing h5-open-106 fixes")