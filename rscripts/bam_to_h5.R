suppressMessages(library(getopt, quietly = T))
# Determine location of provenance.R relative to current file
source(file.path(dirname(getopt::get_Rscript_filename()), "provenance.R"))
suppressMessages(library(Rsamtools, quietly = T))
suppressMessages(library(rtracklayer, quietly = T))
suppressMessages(library(rhdf5, quietly = T))
suppressMessages(library(parallel, quietly = T))
suppressMessages(library(optparse,quietly = T))
suppressMessages(library(RcppRoll, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(tidyr, quietly = T))

reads_to_list <- function(gene_location, bam_file, read_range, left_flank,right_flank, mult_exon=TRUE)
{
  # Computes matrix of read counts with starting position and read length
  # on a signle gene whose co-ordinates are supplied as arguments
  #  
  # Args:
  #   gene_location: (list) gff line containing gene co-ordinates
  #   bam_file: (string) name of bam_file with reads in it
  #   read_range: (int) vector of integer read lengths, e.g. 15:50 
  #   flank: (int) width of flanking region to include outside gene body/CDS
  #   mult_exon: If TRUE, uses only the first supplied exon for gene co-ordinates
  #
  # Returns:
  #   output: matrix of integer counts for each row/start position and column/read length
  #

  if(!mult_exon)
    # If the gene is a single exon gene but multiple GRanges specified,
    # use the first GRange
  { gene_location <- gene_location[1]
  }
  
  # Number of read length categories
  read_range_len <- length(read_range)

  # Specify output matrix to store read data
  gene_length = sum(sum(coverage(gene_location)))
  output <- matrix(0,nrow=read_range_len,ncol=(gene_length + left_flank + right_flank))
  
  # If CDS is on negative strand
  if(as.character(strand(gene_location)[1])=="-"){
    tmp <- left_flank
    left_flank <- right_flank
    right_flank <- tmp
  }
  
  # Expand genomic locations to include flanking region positions
  if(length(gene_location)==1) 
  {	# For single exon genes
    start(gene_location) <- start(gene_location) - left_flank
    end(gene_location) <- end(gene_location) + right_flank
  }else{
    # For multiple exon genes
    start(gene_location)[start(gene_location)==min(start(gene_location))] <- 
        start(gene_location)[start(gene_location)==min(start(gene_location))] - left_flank
    end(gene_location)[end(gene_location)==max(end(gene_location))] <- 
        end(gene_location)[end(gene_location)==max(end(gene_location))] + right_flank
  }
  
  # Specify variables to read in and create a parameter file for bam data
  bam_what <- c("strand", "pos", "qwidth")
  bam_param <- ScanBamParam(which = gene_location, what = bam_what)
  
  # Read in bam data
  bam_data <- scanBam(bam_file, param=bam_param)
  
  # Subset reads that are on the same strand at the genomic location
  read_strand <- unlist(lapply(bam_data,function(x) x$strand ))
  read_location <- unlist(lapply(bam_data,function(x) x$pos ))[read_strand==as.factor(strand(gene_location)[1])]
  read_width <- unlist(lapply(bam_data,function(x) x$qwid ))[read_strand==as.factor(strand(gene_location)[1])]
  
  # Column numbers based on genomic position
  nucleotide_pos <- unlist(which(coverage(gene_location)[seqnames(gene_location)[1]]==1))
  
  # If the specified flanking regions of a gene end up outside the chromosome locations (<0)
  # add pseudo columns with negative numbers
  if(min(start(gene_location))<min(nucleotide_pos))
  { nucleotide_pos <- c(min(start(gene_location)):0,nucleotide_pos)
  }
  
  # Count reads whose 5' ends map to each nucleotide and save them in output matrix
  if(all(strand(gene_location)=="+")){
    # For genes on positive strands
    j <- 1 # counter for output row
    for(i in read_range){
      read_loc_len_i <- read_location[read_width==i] # subset reads of a particular length 
      count_reads_len_i <- table(factor(read_loc_len_i,levels=nucleotide_pos)) # count reads of a particular length 
      output[j,] <- c(count_reads_len_i)
      j <- j+1
    } 
  } 
  
  if(all(strand(gene_location)=="-")){
    # For genes on negative strands
    j <- 1 # counter for output row
    for(i in read_range){
      read_loc_len_i <- read_width[read_width==i]+read_location[read_width==i]-1 # subset reads of a particular length 
      count_reads_len_i <- table(factor(read_loc_len_i,levels=nucleotide_pos)) # count reads of a particular length
      output[j,] <- c(rev(count_reads_len_i))
      j <- j+1
    }
  }
  
  # Return the output matrix
  return(output);
}


# define input options for optparse package
option_list <- list( 
    make_option("--num-processes", type="integer", default=1,
        help="Number of cores for parallelization"),
    make_option("--min-read-length", type="integer", default=10,
        help="Minimum read length in H5 output"),
    make_option("--max-read-length", type="integer", default=50,
        help="Maximum read length in H5 output"),
    make_option("--buffer", type="integer", default=250,
        help="Length of flanking region around the CDS"),
    make_option("--primary-id", type="character", default="gene_id",
        help="Primary gene IDs to access the data (YAL001C, YAL003W, etc.)"),
    make_option("--secondary-id", type="character", default=NULL,
        help="Secondary gene IDs to access the data (COX1, EFB1, etc.)"),
    make_option("--bam-file", type="character", default="input.bam",
        help="Location of BAM input file"),
    make_option("--hd-file", type="character", default="output.h5",
        help="Location of H5 output file"),
    make_option("--orf-gff-file", type="character", default=NULL,
        help="GFF2/GFF3 annotation file"),
    make_option("--dataset", type="character", default="data",
        help="Name of the dataset"),
    make_option("--stop-in-cds", type="logical", default=FALSE,
        help="Are stop codons part of the CDS annotations in GFF?"),
    make_option("--is-riboviz-gff", type="logical", default=TRUE,
                help="Is the GFF file with UTR5, CDS, and UTR3 elements per gene?"),
    make_option("--is-test-run", type="logical", default=FALSE,
        help="Is this a test run?")
    )

print_provenance(get_Rscript_filename())
# Read in commandline arguments
opt <- parse_args(OptionParser(option_list=option_list),
                  convert_hyphens_to_underscores=TRUE)

if (opt$secondary_id == "NULL") {
  # unquote NULL option in secondary_id
  secondary_id <- NULL
}
attach(opt)

print("bam_to_h5.R running with parameters:")
opt

#####
# defunct alternative options for command line input
# comargs=(commandArgs(TRUE))
# 
# if(length(comargs)==0){
#   print("No arguments supplied.")
#   quit(status=1)
# } else if(length(comargs)==1 & endsWith(comargs[[1]],".yaml")){
#   library(yaml,quietly=T)
#   print("Arguments supplied in yaml file.")
#   configopts <- yaml.load_file(comargs[[1]])
#   attach(configopts)
# } else {
#   print("Arguments supplied in command line.")
#   for(i in 1:length(comargs)){
#     eval(parse(text=comargs[[i]]))
#   }
# }
#####

debugVignette <- FALSE
if (debugVignette) {
    # to debug bam_to_h5.R in Riboviz vignette
    # supplies non-default options needed.
    bam_file <- 'vignette/output/SRR1042855_s1mi.bam'
    hd_file <- 'vignette/output/SRR1042855_s1mi.h5'
    orf_gff_file <- 'vignette/input/yeast_YAL_CDS_w_250utrs.gff3'
    num_processes <- 4
    is_test_run <- TRUE
    primary_id <- "Name"
    secondary_id <- NULL
    # gene <- "YAL038W" # use this to debug single gene for reads_to_list
}



# Range of read lengths 
read_range <- min_read_length:max_read_length

# Read in the positions of all exons/genes in GFF format and subset CDS locations
gff <- readGFFAsGRanges(orf_gff_file)

if(!is_riboviz_gff){
  gff <- gff[gff$type=="CDS"]
}

# Read in the list of genes
genes <- unique(mcols(gff)[primary_id][,1])
if(!is.null(secondary_id)){
  alt_genes <- as.list(unique(mcols(gff)[secondary_id][,1]))
  names(alt_genes) <- genes
}
gff_pid <- mcols(gff)[primary_id][,1]

print("Mapping reads to CDS")

if(is_test_run){
  genes <- genes[1:50]
  gff <- gff[ gff_pid %in% genes]
}

## desirable to check seqlevels here rather than wait for later error
# get seqinfo from bam_file
# bam_fileSeqInfo <- seqinfo(scanBamHeader(bam_file))

# Map the reads to individual nucleotide position for each gene
outputList <- 
     mclapply(genes, # parallel version, take out for debugging
    # lapply(genes,
                       function(gene){
                           # select gene location
                           gene_location=gff[gff_pid==gene]
                           # restrict seqlevels to avoid scanBam error
                           geneseqname <- 
                               as.character(seqnames(gene_location)[1])
                           seqlevels(gene_location) <- geneseqname
                           # if GFF contains UTR5, CDS, and UTR3 elements
                           if(is_riboviz_gff){
                             Buffer_l <- width(gene_location[gene_location$type=="UTR5"])
                             Buffer_r <- width(gene_location[gene_location$type=="UTR3"])
                             gene_location <- gene_location[gene_location$type=="CDS"]
                           }else{
                             Buffer_l <- buffer
                             Buffer_r <- buffer
                           }
                           # get reads to list
                         reads_to_list(gene_location=gene_location, 
                                       bam_file=bam_file, 
                                       read_range=read_range, 
                                       left_flank=Buffer_l,
                                       right_flank=Buffer_r, 
                                       mult_exon=TRUE)
                         }
                    # )
 , 
                        mc.cores=num_processes)

names(outputList) <- genes

# HDF5 section
print("Saving mapped reads in a H5 file")

h5createFile(hd_file) # Create the output h5 file
fid <- H5Fopen(hd_file) # Filehandle for the h5 file

# Start codon position
if(!is_riboviz_gff){
  start_cod <- (buffer+1):(buffer+3)
}
# Stop codon offset
if(stop_in_cds){
  offset <- 2
}else{
  offset <- -1
}

# To create symbolic links for alternate gene ids
if(!is.null(secondary_id)){
  base_gid <- H5Gopen(fid,"/")
}

for(gene in genes){
  # Get the output matrix of read counts by position and length for a gene 
  output <- outputList[[gene]]
  
  # Location of start and stop codon nucleotides in output matrix
  gene_location=gff[gff_pid==gene]
  if(is_riboviz_gff){
    start_codon_loc <- start(gene_location[gene_location$type=="CDS"])
    start_cod <- start_codon_loc:(start_codon_loc+2)
    stop_codon_loc <- start(gene_location[gene_location$type=="UTR3"])-3
  }else{
    stop_codon_loc <- ncol(output)-buffer-offset
  }
  
  stop_cod <- stop_codon_loc:(stop_codon_loc+2)
  
  # Create H5 groups for each gene
  h5createGroup(fid,gene)
  h5createGroup(fid,paste(gene,dataset,sep="/"))
  h5createGroup(fid,paste(gene,dataset,"reads",sep="/"))

  mapped_reads <- paste(gene,dataset,"reads",sep="/")
  
  # Symbolic link with alternate ids
  if(!is.null(secondary_id)){ 
    if(alt_genes[[gene]]!=gene){
      H5Lcreate_external(hd_file, gene, base_gid, alt_genes[[gene]])
    }
  }
  
  # Create a handle for H5 group for the gene
  gid <- H5Gopen(fid, mapped_reads)
  
  # Specify attributes of the gene
  h5createAttribute(gid,"reads_total",c(1,1))
  h5createAttribute(gid,"buffer_left",c(1,1))
  h5createAttribute(gid,"buffer_right",c(1,1))
  h5createAttribute(gid,"start_codon_pos",c(1,3))
  h5createAttribute(gid,"stop_codon_pos",c(1,3))
  h5createAttribute(gid,"reads_by_len",c(1,length(read_range)))
  h5createAttribute(gid,"lengths",c(1,length(read_range)))
  
  # # Count number of reads that map to the CDS of a gene
  # # We include reads mapping to 25 bases up/downstream of CDS in this count
  # ifelse(start_cod[1]>25, up_atg <- (start_cod[1]-25), up_atg <- start_cod[1])
  # ifelse(ncol(output)>(stop_cod[3]+25), down_taa <- (stop_cod[3]+25), down_taa <- stop_cod[3])
  # 
  # cod_total <- output[,up_atg:down_taa]
  
  # Write the attributes of the gene group
  h5writeAttribute.integer(sum(output),gid,name="reads_total")
  h5writeAttribute.integer((start_codon_loc-1),gid,name="buffer_left")
  h5writeAttribute.integer((ncol(output)-stop_cod[3]),gid,name="buffer_right")
  h5writeAttribute.integer(start_cod,gid,name="start_codon_pos")
  h5writeAttribute.integer(stop_cod,gid,name="stop_codon_pos")
  h5writeAttribute.integer(read_range,gid,name="lengths")
  h5writeAttribute.integer(apply(output,1,sum),gid,name="reads_by_len")
  
  # Specify a dataset within the gene group to store the values and degree of compression
  read_data <- paste(mapped_reads,"data",sep="/")
  h5createDataset(fid,read_data,dim(output), storage.mode="integer", chunk=c(1,ncol(output)),level=7)
  
  # Write the dataset within the gene group
  h5write(output,fid,name=read_data,start=c(1,1))
  
  # Close H5 gene group handle
  H5Gclose(gid)
}

if(!is.null(secondary_id)){
  H5Gclose(base_gid)
}

H5close()

print("bam_to_h5.R done")
