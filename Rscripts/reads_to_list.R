reads_to_list <- function(gene, gene_location, bamFile, read_range=10:50, left_buffer=250, right_buffer=247,mult_exon=TRUE)
{
  flank <- right_buffer+3
  
  if(!mult_exon)
  { gene_location <- gene_location[1]
  }
	read_range_len <- length(read_range)
	start_cod <- (left_buffer+1):(left_buffer+3)

	# Specify output matrix
	output <- matrix(0,nrow=read_range_len,ncol=(sum(sum(coverage(gene_location)))+2*flank))

	# Check for introns
	if(length(gene_location)==1)
	{	start(gene_location) <- start(gene_location)-flank
		end(gene_location) <- end(gene_location)+flank
	}else{
		start(gene_location)[start(gene_location)==min(start(gene_location))] <- start(gene_location)[start(gene_location)==min(start(gene_location))]-flank
		end(gene_location)[end(gene_location)==max(end(gene_location))] <- end(gene_location)[end(gene_location)==max(end(gene_location))]+flank
	}

	# Read in bam data
	what <- c("strand", "pos", "qwidth")
	param <- ScanBamParam(which = gene_location, what = what)
	bam <- scanBam(bamFile, param=param)

	read_strand <- unlist(lapply(bam,function(x)x$strand))
	read_location <- unlist(lapply(bam,function(x)x$pos))[read_strand==as.factor(strand(gene_location)[1])]
	read_width <- unlist(lapply(bam,function(x)x$qwid))[read_strand==as.factor(strand(gene_location)[1])]

	# Column numbers based on genomic position
	column_pos <- unlist(which(coverage(gene_location)[seqnames(gene_location)[1]]==1))
	if(start(gene_location)<min(column_pos))
  { column_pos <- c(start(gene_location):0,column_pos)
	}

	if(all(strand(gene_location)=="+"))
	{	j <- 1
		for(i in read_range)
		{	x <- read_location[read_width==i]
			ty <- table(factor(x,levels=column_pos))
			output[j,] <- c(ty)
			j <- j+1
		} 
	} 
	if(all(strand(gene_location)=="-"))
	{	j <- 1
		for(i in read_range)
		{	x <- read_width[read_width==i]+read_location[read_width==i]-1
			ty <- table(factor(x,levels=column_pos))
			output[j,] <- c(rev(ty))
			j <- j+1
		}
	}

	return(output);
}
