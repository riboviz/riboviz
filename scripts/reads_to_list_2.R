reads_to_list <- function(gene_location, bamFile, read_range, flank, mult_exon=TRUE)
{
  # If the gene is a single exon gene but multiple GRanges specified,
  # use the first GRange
  if(!mult_exon)
  { gene_location <- gene_location[1]
  }
  
  # Number of read length categories
  read_range_len <- length(read_range)

  # Specify output matrix to store read data
  output <- matrix(0,nrow=read_range_len,ncol=(sum(sum(coverage(gene_location)))+2*flank))
  
  # Expand genomic locations to include flanking region positions
  if(length(gene_location)==1) 
  {	# For single exon genes
    start(gene_location) <- start(gene_location)-flank
    end(gene_location) <- end(gene_location)+flank
  }else{
    # For multiple exon genes
    start(gene_location)[start(gene_location)==min(start(gene_location))] <- start(gene_location)[start(gene_location)==min(start(gene_location))]-flank
    end(gene_location)[end(gene_location)==max(end(gene_location))] <- end(gene_location)[end(gene_location)==max(end(gene_location))]+flank
  }
  
  # Specify variables to read in and create a parameter file for bam data
  what <- c("strand", "pos", "qwidth")
  param <- ScanBamParam(which = gene_location, what = what)
  
  # Read in bam data
  bam_data <- scanBam(bamFile, param=param)
  
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
  # For genes on positive strands
  if(all(strand(gene_location)=="+")){
    j <- 1 # counter for output row
    for(i in read_range){
      read_loc_len_i <- read_location[read_width==i] # subset reads of a particular length 
      count_reads_len_i <- table(factor(read_loc_len_i,levels=nucleotide_pos)) # count reads of a particular length 
      output[j,] <- c(count_reads_len_i)
      j <- j+1
    } 
  } 
  
  # For genes on negative strands
  if(all(strand(gene_location)=="-")){
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
