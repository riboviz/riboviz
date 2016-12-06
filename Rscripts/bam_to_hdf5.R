# Read in dependent packages
library(Rsamtools, quietly=T)
library(rtracklayer, quietly=T)
library(rhdf5, quietly=T)

# Initialize parameters
ncores=NULL
read_range <- 15:50
left_buffer <- 250 
right_buffer <- 247 

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

# Read in mapping function
source("/home/prshah/codes/ribosome_profiling/R/reads_to_list.R")

# Read in the positions of all exons/genes in GTF/GFF format
gtf <- import(con=gtfFile, format="gtf")

# Read in the list of genes
genes <- unique(gtf$Name[gtf$type=="CDS"])

# Map the reads to individual nucleotide position for each gene
if(is.null(ncores))
{ outputList <- lapply(genes, function(x){gene_location <- gtf[gtf$Name==x & gtf$type=="CDS"];reads_to_list(gene=x, gene_location=gene_location, bamFile=bamFile, mult_exon=TRUE);})
}else{
  library(parallel, quietly=T)
  outputList <- mclapply(genes, function(x){gene_location <- gtf[gtf$Name==x & gtf$type=="CDS"];reads_to_list(gene=x, gene_location=gene_location, bamFile=bamFile, mult_exon=TRUE);}, mc.cores=ncores)
}
names(outputList) <- genes

# HDF5 section

h5createFile(hdfile) # Create an .h5 file
fid <- H5Fopen(hdfile) # Filehandle for the .h5 file

read_range_len <- length(read_range)
start_cod <- (left_buffer+1):(left_buffer+3)

for(gene in genes)
{
	output <- outputList[[gene]]
	
	# HDF5 section
	mapped_reads <- paste(gene,dataset,"reads",sep="/")

	h5createGroup(fid,gene)
	h5createGroup(fid,paste(gene,dataset,sep="/"))
	h5createGroup(fid,paste(gene,dataset,"reads",sep="/"))

	x <- ncol(output)-right_buffer-2
	stop_cod <- x:(x+2)
	
	gid <- H5Gopen(fid, mapped_reads)
	
	h5createAttribute(gid,"reads_total",c(1,1))
	h5createAttribute(gid,"buffer_left",c(1,1))
	h5createAttribute(gid,"buffer_right",c(1,1))
	h5createAttribute(gid,"start_codon_pos",c(1,3))
	h5createAttribute(gid,"stop_codon_pos",c(1,3))
	h5createAttribute(gid,"reads_by_len",c(1,length(read_range)))
	h5createAttribute(gid,"lengths",c(1,length(read_range)))
	
	cod_total <- output[,(start_cod[1]-25):(stop_cod[3]+25)]

	h5writeAttribute.integer(sum(cod_total),gid,name="reads_total")
	h5writeAttribute.integer(left_buffer,gid,name="buffer_left")
	h5writeAttribute.integer(right_buffer,gid,name="buffer_right")
	h5writeAttribute.integer(start_cod,gid,name="start_codon_pos")
	h5writeAttribute.integer(stop_cod,gid,name="stop_codon_pos")
	h5writeAttribute.integer(read_range,gid,name="lengths")
	h5writeAttribute.integer(apply(cod_total,1,sum),gid,name="reads_by_len")
	
	read_data <- paste(mapped_reads,"data",sep="/")
	h5createDataset(fid,read_data,dim(output), storage.mode="integer", chunk=c(1,ncol(output)),level=7)
	h5write(output,fid,name=read_data,start=c(1,1))
	
	H5Gclose(gid)
}
H5close()
