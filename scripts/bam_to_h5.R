# Read in dependent packages
library(Rsamtools, quietly=T)
library(rtracklayer, quietly=T)
library(rhdf5, quietly=T)
library(parallel, quietly=T)

reads_to_list <- function(gene_location, bamFile, read_range, flank, mult_exon=TRUE)
{
  # Computes matrix of read counts with starting position and read length
  # on a signle gene whose co-ordinates are supplied as arguments
  #  
  # Args:
  #   gene_location: (list) gff line containing gene co-ordinates
  #   bamFile: (string) name of bamFile with reads in it
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
  gene_length = sum(sum(coverage(gene_location))
  output <- matrix(0,nrow=read_range_len,ncol=gene_length + 2*flank))
  
  # Expand genomic locations to include flanking region positions
  if(length(gene_location)==1) 
  {	# For single exon genes
    start(gene_location) <- start(gene_location) - flank
    end(gene_location) <- end(gene_location) + flank
  }else{
    # For multiple exon genes
    start(gene_location)[start(gene_location)==min(start(gene_location))] <- 
        start(gene_location)[start(gene_location)==min(start(gene_location))] - flank
    end(gene_location)[end(gene_location)==max(end(gene_location))] <- 
        end(gene_location)[end(gene_location)==max(end(gene_location))] + flank
  }
  
  # Specify variables to read in and create a parameter file for bam data
  bam_what <- c("strand", "pos", "qwidth")
  bam_param <- ScanBamParam(which = gene_location, what = bam_what)
  
  # Read in bam data
  bam_data <- scanBam(bamFile, param=bam_param)
  
  # Subset reads that are on the same strand at the genomic location
  read_strand <- unlist(lapply(bam_data,function(x)x$strand))
  read_location <- unlist(lapply(bam_data,function(x)x$pos))[read_strand==as.factor(strand(gene_location)[1])]
  read_width <- unlist(lapply(bam_data,function(x)x$qwid))[read_strand==as.factor(strand(gene_location)[1])]
  
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


# Initialize default parameters
Ncores <- 1 # Number of cores for parallelization
MinReadLen <- 10 # Minimum read length in H5 output
MaxReadLen <- 50 # Maximum read length in H5 output
Buffer <- 250 # Length of flanking region around the CDS
PrimaryID <- "gene_id" # Primary gene IDs to access the data (YAL001C, YAL003W, etc.)
SecondID <- NULL # Secondary gene IDs to access the data (COX1, EFB1, etc.)
bamFile <- "input.bam" # Location of BAM file
hdfile <- "output.h5" # Location of H5 output file
gffFile <- NULL # Location of GFF2/GFF3 annotation file
dataset <- "data" # Name of the dataset
StopInCDS <- FALSE # Are stop codons part of the CDS annotations in GFF?
test <- FALSE # Is this a test run

# Read in commandline arguments
args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
  quit(status=1)
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

# Range of read lengths 
read_range <- MinReadLen:MaxReadLen

# Read in the positions of all exons/genes in GFF format and subset CDS locations
gff <- readGFFAsGRanges(gffFile)
gff <- gff[gff$type=="CDS"]

# Read in the list of genes
genes <- unique(mcols(gff)[PrimaryID][,1])
if(!is.null(SecondID)){
  alt_genes <- as.list(unique(mcols(gff)[SecondID][,1]))
  names(alt_genes) <- genes
}
gff_pid <- mcols(gff)[PrimaryID][,1]

print("Mapping reads to CDS")

if(test){
  genes <- genes[1:50]
}

# Map the reads to individual nucleotide position for each gene
outputList <- mclapply(genes, 
                       function(x){
                         gene_location <- gff[gff_pid==x]
                         reads_to_list(gene_location=gene_location, 
                                       bamFile=bamFile, 
                                       read_range=read_range, 
                                       flank=Buffer, 
                                       mult_exon=TRUE);
                         }, 
                       mc.cores=Ncores)

names(outputList) <- genes

# HDF5 section
print("Saving mapped reads in a H5 file")

h5createFile(hdfile) # Create the output h5 file
fid <- H5Fopen(hdfile) # Filehandle for the h5 file

# Start codon position
start_cod <- (Buffer+1):(Buffer+3)
# Stop codon offset
if(StopInCDS){
  offset <- 2
}else{
  offset <- -1
}

# To create symbolic links for alternate gene ids
if(!is.null(SecondID)){
  base_gid <- H5Gopen(fid,"/")
}

for(gene in genes){
  # Get the output matrix of read counts by position and length for a gene 
  output <- outputList[[gene]]
  
  # Create H5 groups for each gene
  h5createGroup(fid,gene)
  h5createGroup(fid,paste(gene,dataset,sep="/"))
  h5createGroup(fid,paste(gene,dataset,"reads",sep="/"))

  mapped_reads <- paste(gene,dataset,"reads",sep="/")
  
  # Symbolic link with alternate ids
  if(!is.null(SecondID) & alt_genes[[gene]]!=gene){
    H5Lcreate_external(hdfile, gene, base_gid, alt_genes[[gene]])
  }
  
  # Location of stop codon nucleotides in output matrix
  stop_codon_loc <- ncol(output)-Buffer-offset
  stop_cod <- stop_codon_loc:(stop_codon_loc+2)
  
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
  
  # Count number of reads that map to the CDS of a gene
  # We include reads mapping to 25 bases up/downstream of CDS in this count
  cod_total <- output[,(start_cod[1]-25):(stop_cod[3]+25)]
  
  # Write the attributes of the gene group
  h5writeAttribute.integer(sum(cod_total),gid,name="reads_total")
  h5writeAttribute.integer(Buffer,gid,name="buffer_left")
  h5writeAttribute.integer(Buffer,gid,name="buffer_right")
  h5writeAttribute.integer(start_cod,gid,name="start_codon_pos")
  h5writeAttribute.integer(stop_cod,gid,name="stop_codon_pos")
  h5writeAttribute.integer(read_range,gid,name="lengths")
  h5writeAttribute.integer(apply(cod_total,1,sum),gid,name="reads_by_len")
  
  # Specify a dataset within the gene group to store the values and degree of compression
  read_data <- paste(mapped_reads,"data",sep="/")
  h5createDataset(fid,read_data,dim(output), storage.mode="integer", chunk=c(1,ncol(output)),level=7)
  
  # Write the dataset within the gene group
  h5write(output,fid,name=read_data,start=c(1,1))
  
  # Close H5 gene group handle
  H5Gclose(gid)
}
if(!is.null(SecondID)){
  H5Gclose(base_gid)
}

H5close()
