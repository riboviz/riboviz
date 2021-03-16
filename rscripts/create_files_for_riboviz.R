## Authors: Alexander Cope, Premal Shah
## An R script for generating riboviz-style CDS and GFF3 files from a genome FASTA (i.e. sequences represent chromosomes, plasmids, etc. but not transcripts, CDS, etc.) 
## and the corresponding GFF3 file. Can also generate other optional files (H5 files, codon position .Rdata files) for use with the
## riboviz pipeline. Note that this script assumes sequence names in FASTA and GFF3 file (left-most column) match. This may require some 
## independent pre-processing of the files by the user. 

library(argparse)
library(Biostrings)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(parallel)
library(rhdf5)

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="Input DNA sequences. Should contain genome (i.e. sequence of each chromosome) or transcripts. Should be file path",type="character")
parser$add_argument("-g","--gff",help="GFF3 file corresponding to input DNA sequences.",type="character")
parser$add_argument("--out_dir",help="Output directory for Riboviz-style CDS and GFF files. Will be created if does not exist.",type="character",default="./")
parser$add_argument("--out_cds",help="Name of Riboviz-style CDS file",type="character",default="riboviz_cds.fasta")
parser$add_argument("--out_gff",help="Name of Riboviz-style GFF3 file",type="character",default="riboviz_gff.gff3")
parser$add_argument("-s","--seq_id",help="GFF column to use as the sequence ids",type="character",default="Name")
parser$add_argument("-b","--buffer",help="Buffer to use for UTRs",type="integer",default=250)
parser$add_argument("--h5_file",help="File name for createing H5 file. If not initialized, file will not be created.",type="character",default=NULL)
parser$add_argument("--codon_data_file",help="File name for codon position .Rdata file. If not initialized, file will not be created.",type="character",default=NULL)
parser$add_argument("--num_cores",help="Number of cores to use for parallelizable processes.",type="integer",default=1)
parser$add_argument("--codons_exclude",help="Exclude the first n codons when creating codon_data_file, where n is specified by this argument",default=0)
parser$add_argument("--remove_trailing",help="Remove trailing info from names to be used for CDS, e.g. remove anything after '_' or '.'",type="character",default="_|\\.")
parser$add_argument("--filter_seq",help="A comma-separated list of filtering criteria to apply to the GFF3 file, e.g. 'type:CDS,orf_classification:Verified,orf_classification:Uncharacterized'. Use 'notNA' to filter values that are NA, e.g. 'orf_classification:!NA'.",type="character",default="type:CDS")
parser$add_argument("--exons_preordered",help="Some GFF3 files have exons pre-ordered such that exon with start codon is listed first. Effects how multi-exon genes will be combined.",action="store_true")

args <- parser$parse_args()
input <- args$input
gff <- args$gff
output_dir <- args$out_dir
output_cds <- args$out_cds
output_gff <- args$out_gff
seq_id <- args$seq_id
buffer <- args$buffer
h5_file <- args$h5_file
codon_data_file <- args$codon_data_file
num_cores <- args$num_cores
codons_exclude <- args$codons_exclude
remove_trailing <- args$remove_trailing
filter_seq <- args$filter_seq
exons_preordered <- args$exons_preordered

##### Helper functions #########################################################




# Funciton to resize GRanges while keeping track of multi-exon genes
resize_gff_flank <- function(gene_location,flank){
  if(length(gene_location)==1){	
    start(gene_location) <- start(gene_location)-flank
    end(gene_location) <- end(gene_location)+flank
  }else{
    start(gene_location)[start(gene_location)==min(start(gene_location))] <- start(gene_location)[start(gene_location)==min(start(gene_location))]-flank
    end(gene_location)[end(gene_location)==max(end(gene_location))] <- end(gene_location)[end(gene_location)==max(end(gene_location))]+flank
  }
  return(gene_location)
}

#' A function Premal used for getting the UTRs by comparing the exon and CDS.
#' This seemed to be based on some deprecated functionality in IRanges, hence
#' the use of IRanges:::unlist_as_integer
#' I think a better strategy going forward is to just use the five_prime_UTR
#' and three_prime_UTR annotations.
get_utr_cds_len <- function(exon,cds,strand){
  posdiff <- IRanges:::unlist_as_integer(setdiff(ranges(exon),ranges(cds)))
  utr5 <- sum(posdiff<min(start(cds)))
  clen <- sum(width(cds))
  utr3 <- sum(posdiff>max(end(cds)))
  if(strand=="-"){
    tmp <- utr3
    utr3 <- utr5
    utr5 <- tmp
  }
  return(c(utr5,clen,utr3))
}


# Function to find positions for a codon in a CDS
codon_pos_locator <- function(seq,codon,gene){
  out <- NULL
  pos <- which(seq==codon)
  if(length(pos)>=1)
    out <- cbind(gene,pos)
  return(out)
}

##### End Helper Functions #####################################################

#' A function for filtering the GFF file based on some criteria
#' Primarily, we expect to filter CDS files, but other situations might arise.
#' For example, SGD has orf_classification as "Verified", "Uncharacterized", 
#' and "Dubious". 
#' @param gff A GRanges object from a GFF3 file
#' @param criteria A comma-separated string representing the criteria to filter of the form <column:criteria>, 
#'                 e.g. type:CDS, orf_classification:Verified.
#' @return Updated GRanges object
filterGFFByCriteria <- function(gff,criteria)
{
  criteria <- strsplit(unlist(strsplit(criteria,",")),":")
  df_criteria <- as.data.frame(do.call(rbind,criteria))
  criteria_keys <- unique(df_criteria[,1])
  for (i in criteria_keys)
  {
    criteria_conditions <- df_criteria[which(df_criteria[,1] == i),2]
    criteria_conditions_not_na <- criteria_conditions[which(criteria_conditions != "notNA")]
    gff <- gff[mcols(gff)[[i]] %in% criteria_conditions_not_na]
    ## If want to also remove NA values, need to do this using is.na
    if (length(criteria_conditions_not_na) < length(criteria_conditions))
    {
      gff <- gff[!is.na(mcols(gff)[[i]])]
    }
  }
  return(gff)
}

#' A function to remove trailing information from gene names, such as version ids
#' @param gff A GRanges object from a GFF3 file
#' @param column The column from which to remove trailing ids
#' @param regex A regular expression to trim from the ids
#' @return Updated GRanges object
removeTrailing <- function(gff,column,regex="_|\\.|\\:")
{
  gff$Name <- sapply(strsplit(mcols(gff)[[column]],regex), `[`, 1)
  return(gff)
}

#' Read in a genome file (not a transcriptome file).
#' @param genome.file File path leading to 
#' @return Return genome as DNAStringSet object
readInGenomeFasta <- function(genome.file)
{
  genome <- readDNAStringSet(genome.file,format = "fasta")
  return(genome)
}


#' Given a buffer to use for fixed UTRs, append to ends of chromosomes and shift GRanges obhect
#' @param genome DNAStringSet object with genome
#' @param annot GRanges object with GFF3
#' @param buffer Size for UTRs
#' @return Return list with genome (as DNAStringSet object) and annotations/GFF3 (as GRanges object)
addBufferToGenome <- function(genome,annot,buffer=250)
{
 genome <- xscat(DNAString(paste(rep("N",buffer),collapse="")),
                        genome,
                        DNAString(paste(rep("N",buffer),collapse="")))
 annot <- shift(annot,buffer)
 return(list(Genome=genome,Annotation=annot))
}


#' Given a genome and a set of names, return the genome with the new names. Note that the assumes new.names is in appropriate order
#' @param genome DNAStringSet object with genome
#' @param new.names New names to give genomes
#' @return Updated genome as DNAStringSet object
renameGenome <- function(genome,new.names)
{
  if (length(new.names) != length(genome))
  {
    print("The number of names provided does not match the length of the genome object. Genome will be returned unmodified.") 
  } else{
    names(genome) <- new.names
  }
  return(genome)
}


#' Given a genome, a GFF annotation, and a buffer, create A riboviz-style CDS.
#' @param genome DNAStringSet object with genome
#' @param gff A GRanges object from a GFF3 file
#' @param buffer Size of buffers for UTRs
#' @param exons_preordered Do exons need to be reordered? Some GFF3 files have exons for multi-exon genes pre-ordered such that the exon with the start codon is first in the file.
#' @return CDS as DNAStringSet object
convertGenomeToCDSFile <- function(genome,gff_annot,buffer,exons_preordered = FALSE)
{
  tmp_gff <- sapply(unique(gff_annot$Name),
                   function(x){
                     gene_location <- gff_annot[gff_annot$Name==x]; 
                     resize_gff_flank(gene_location = gene_location, flank=buffer)
                   })
  cds_flank_annot <- unlist(GRangesList(tmp_gff))
  
  # Obtain sequences of CDS+UTRs from genome
  cds_flank_seq <- genome[cds_flank_annot]
  names(cds_flank_seq) <- cds_flank_annot$Name
  
  # Combine multi-exon genes into a single transcript
  output_seqlist <- list()
  cc <- 1
  for(i in unique(names(cds_flank_seq))){
    if (!exons_preordered) {
      output_seqlist[[cc]] <- c(unlist(cds_flank_seq[names(cds_flank_seq)==i]))
      if(any(strand(cds_flank_annot[cds_flank_annot$Name==i])=="-")){
        ## In the case of exons not preordered, can just apply reverseComplement to entire sequence. 
        output_seqlist[[cc]] <- reverseComplement(output_seqlist[[cc]])
      }
    } else{
      tmp.seq.list <- cds_flank_seq[names(cds_flank_seq)==i]
      if (any(strand(cds_flank_annot[cds_flank_annot$Name==i])=="-")){
        ## In the case of preordered exons, apply reverseComplement to all exons individually and then merge
        output_seqlist[[cc]] <- c(unlist(reverseComplement(tmp.seq.list)))
      } else {
        output_seqlist[[cc]] <- c(unlist(tmp.seq.list))
      }
    }
    cc <- cc+1
  }
  output_seq <- DNAStringSet(output_seqlist)
  names(output_seq) <- unique(names(cds_flank_seq))
  return(output_seq)
}

#' Given CDS sequences and buffer, generate a riboviz-style GFF
#' @param cds DNAStringSet object with cds sequences 
#' @param gff A GRanges object from a GFF3 file
#' @param buffer Size of buffers for UTRs
#' @return Riboviz-style GFF as GRanges object
createRibovizGFF <- function(cdna,buffer)
{
  len_gff <- c(matrix(c(rep(buffer,length(cdna)), (width(cdna)-(buffer*2)), rep(buffer,length(cdna))),nrow=3,byrow = T))
  start_gff <- c(matrix(c(rep(1,length(cdna)),rep((buffer+1),length(cdna)), (width(cdna)-(buffer-1))),nrow=3,byrow = T))
  type_gff <- rep(c("UTR5","CDS","UTR3"),length(cdna))
  
  # Create a GRange object to be saved as a GFF
  output_gff <- GRanges(seqnames = rep(names(cdna),each=3), 
                        ranges = IRanges(start=start_gff, width=len_gff),
                        strand = "+",type=type_gff,
                        Name=rep(names(cdna),each=3))
  return(output_gff)
}


#' Write riboviz-style CDS
#' @param cds DNAStringSet object with cds sequences 
#' @param output_dir Target directory for outputting the riboviz-style CDS file
#' @param output_cds Name of output file
writeRibovizStyleCDS <- function(cds,output_dir,output_cds)
{
  writeXStringSet(cds,filepath = file.path(output_dir,output_cds),format = "fasta")
}

#' Write riboviz-style GFF3
#' @param gff GRanges object 
#' @param output_dir Target directory for outputting the riboviz-style CDS file
#' @param output_gff Name of output file
writeRibovizStyleGFF <- function(gff,output_dir,output_gff)
{
  export.gff3(riboviz_gff, con=file.path(output_dir,output_gff))
}


#' Creata .Rdata object storing codon positions
#' @param seq CDS sequences as DNAStringSet object
#' @param gff GRanges object 
#' @param codon_position_object File name for .Rdata object storing codon positions
#' @param output_dir Target directory for outputting .Rdata object
#' @param start_pos Codon position data will be stored starting at the nth codon, where n is specified by this parameter. Any sequences less than n will be removed.
#' @param num_cores Specify number of cores to use here
createCodonPositionRData <- function(seq,gff,codon_position_object,output_dir,start_pos=0,num_cores=1)
{
  seq <- seq[gff[gff$type=="CDS"]]
  seq <- seq[width(seq)>start_pos] 
  
  seq <- DNAStringSet(seq,start=start_pos+1) # Trim the first 200 codons from each CDS
  seq <- seq[width(seq)%%3==0] # Ignore any transcripts with frame-shifts
  cods <- sapply(sapply(seq,codons),as.character) # Split the sequences into codons
  
  
  codon_pos <- sapply(names(GENETIC_CODE),
                      function(codon){
                        allpos <- mclapply(names(cods),
                                           function(gene){
                                             codon_pos_locator(seq=cods[[gene]],
                                                               codon=codon,
                                                               gene=gene)
                                           },mc.cores=num_cores)
                        do.call(rbind,allpos)
                      })
  
  codon_pos <- codon_pos[!names(codon_pos) %in% c("TAG","TAA","TGA")]
  # Save output as an RData object
  save(codon_pos, file = file.path(output_dir,codon_position_object))
}


#' Creata .Rdata object storing codon positions
#' @param seq CDS sequences as DNAStringSet object
#' @param gff GRanges object 
#' @param codon.position.object File name for .Rdata object storing codon positions
#' @param output_dir Target directory for outputting .Rdata object
#' @param start Codon position data will be stored starting at the nth codon, where n is specified by this parameter. Any sequences less than n will be removed.
#' @param mc.cores Specify number of cores to use here
createH5File <- function(seq,gff,hdfile,output_dir)
{

  hd.file.path <- file.path(output_dir,hdfile)
  
  seq <- seq[gff[gff$type=="CDS"]] # Restrict sequences to only CDS
  seq <- seq[width(seq)%%3==0] # Ignore any transcripts with frame-shifts
  
  nt_seq <- strsplit(as.character(seq),"")
  cod_seq <- sapply(sapply(seq,codons),as.character) # Split the sequences into codons
  
  # Save seq data as H5 file
  h5createFile(hd.file.path) # Create the output h5 file
  fid <- H5Fopen(hd.file.path) # Filehandle for the h5 file
  
  for(gene in names(nt_seq)){
    # Get the output matrix of read counts by position and length for a gene 
    tmp_nt <- nt_seq[[gene]]
    tmp_cod <- cod_seq[[gene]]
    
    # Create H5 groups for each gene
    h5createGroup(fid,gene)
    
    # Specify a dataset within the gene group to store the values and degree of compression
    nt_dname <- paste(gene,"nt",sep="/")
    cod_dname <- paste(gene,"codon",sep="/")
    
    h5createDataset(fid,nt_dname,dims=c(1,length(tmp_nt)), storage.mode="character",size=2,level=7)
    h5createDataset(fid,cod_dname,dims=c(1,length(tmp_cod)), storage.mode="character",size=4,level=7)
    
    # Write the dataset within the gene group
    h5write(as.matrix(t(tmp_nt)),fid,name=nt_dname,start=c(1,1))
    h5write(as.matrix(t(tmp_cod)),fid,name=cod_dname,start=c(1,1))
  }
  H5close()
}


if (!file.exists(input)){
  stop("FASTA file specified by --input argument not found.")
}
if (!file.exists(gff)){
  stop("GFF file specified by --gff argument not found.")
}

if(!dir.exists(output_dir)){
  print("Output directory not found. Creating...")
  dir.create(output_dir)
}

print("Reading in Genome...")
genome <- readInGenomeFasta(input)
print("Done")
print("Reading in GFF3...")
annot <- readGFFAsGRanges(gff)
print("Done")
current.names <- names(genome)

print("Adding Buffer to genome...")
genome_annot <- addBufferToGenome(genome,annot,buffer=buffer)
print("Done")
genome <- genome_annot$Genome
annot <- genome_annot$Annotation
names(genome) <- current.names

print("Filtering GFF3 by provided criteria...")
annot <- filterGFFByCriteria(gff = annot,criteria = filter_seq)
print("Done")
if (!is.null(remove_trailing))
{
  print("Removing trailing values on gene ids...")
  annot <- removeTrailing(gff = annot, column=seq_id, regex = remove_trailing)
  print("Done")
}
print("Creating riboviz-style CDS...")
gff_annot <- annot

#cds <- convertGenomeToCDSFile(genome,annot,buffer,exons_preordered)
tmp_gff <- sapply(unique(gff_annot$Name),
                  function(x){
                    gene_location <- gff_annot[gff_annot$Name==x]; 
                    resize_gff_flank(gene_location = gene_location, flank=buffer)
                  })
cds_flank_annot <- unlist(GRangesList(tmp_gff))

# Obtain sequences of CDS+UTRs from genome
cds_flank_seq <- genome[cds_flank_annot]
names(cds_flank_seq) <- cds_flank_annot$Name

# Combine multi-exon genes into a single transcript
output_seqlist <- list()
cc <- 1
for(i in unique(names(cds_flank_seq))){
  if (!exons_preordered) {
    output_seqlist[[cc]] <- c(unlist(cds_flank_seq[names(cds_flank_seq)==i]))
    if(any(strand(cds_flank_annot[cds_flank_annot$Name==i])=="-")){
      ## In the case of exons not preordered, can just apply reverseComplement to entire sequence. 
      output_seqlist[[cc]] <- reverseComplement(output_seqlist[[cc]])
    }
  } else{
    tmp.seq.list <- cds_flank_seq[names(cds_flank_seq)==i]
    if (any(strand(cds_flank_annot[cds_flank_annot$Name==i])=="-")){
      ## In the case of preordered exons, apply reverseComplement to all exons individually and then merge
      output_seqlist[[cc]] <- c(unlist(reverseComplement(tmp.seq.list)))
    } else {
      output_seqlist[[cc]] <- c(unlist(tmp.seq.list))
    }
  }
  cc <- cc+1
}
output_seq <- DNAStringSet(output_seqlist)
names(output_seq) <- unique(names(cds_flank_seq))
cds <- output_seq
print("Done")
print("Writing riboviz-style CDS to file...")
writeRibovizStyleCDS(cds,output_dir,output_cds)
print("Done")

print("Creating riboviz-style GFF3...")
riboviz_gff <- createRibovizGFF(cds,buffer = buffer)
print("Done")
print("Writing riboviz-style GFF3 to file...")
writeRibovizStyleGFF(riboviz_gff,output_dir,output_gff)
print("Done")

if (!is.null(codon_data_file))
{
  print("Creating codon position .Rdata file...")
  createCodonPositionRData(seq = cds, gff = riboviz_gff, codon_position_object = codon_data_file, output_dir = output_dir,start_pos = codons_exclude, num_cores = num_cores)
  print("Done")
}

if (!is.null(h5_file))
{
  print("Creating H5 file...")
  createH5File(seq = cds,gff = riboviz_gff, hdfile = h5_file,output_dir=output_dir)
  print("Done")
}
