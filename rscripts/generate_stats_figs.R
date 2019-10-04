# read in dependent packages
suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))
suppressMessages(library(rhdf5))
suppressMessages(library(parallel))
suppressMessages(library(optparse))
suppressMessages(library(RcppRoll))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))

# set ggplot2 theme for plots drawn after this; use dark on light theme
ggplot2::theme_set(theme_bw())

# define input options for optparse package
option_list <- list(
  make_option("--Ncores",
    type = "integer", default = 1,
    help = "Number of cores for parallelization"
  ),
  make_option("--MinReadLen",
    type = "integer", default = 10,
    help = "Minimum read length in H5 output"
  ),
  make_option("--MaxReadLen",
    type = "integer", default = 50,
    help = "Maximum read length in H5 output"
  ),
  make_option("--Buffer",
    type = "integer", default = 250,
    help = "Length of flanking region around the CDS"
  ),
  make_option("--PrimaryID",
    type = "character", default = "gene_id",
    help = "Primary gene IDs to access the data (YAL001C, YAL003W, etc.)"
  ),
  make_option("--hdFile",
    type = "character", default = "output.h5",
    help = "Location of H5 output file"
  ),
  make_option("--dataset",
    type = "character", default = "data",
    help = "Name of the dataset"
  ),
  make_option("--out_prefix",
    type = "character", default = "out",
    help = "Prefix for output files"
  ),
  make_option("--orf_fasta",
    type = "character", default = FALSE,
    help = "FASTA file with nt seq"
  ),
  make_option("--rpf",
    type = "logical", default = TRUE,
    help = "Is the dataset an RPF or mRNA dataset?"
  ),
  make_option("--dir_out",
    type = "character", default = "./",
    help = "Output directory"
  ),
  make_option("--t_rna",
    type = "character", default = NA,
    help = "tRNA estimates in .tsv file"
  ),
  make_option("--codon_pos",
    type = "character", default = NA,
    help = "Codon positions in each gene in .Rdata file"
  ),
  make_option("--orf_gff_file",
    type = "character", default = NA,
    help = "riboviz generated GFF2/GFF3 annotation file"
  ),
  make_option("--features_file",
    type = "character", default = NA,
    help = "features file, columns are gene features and rows are genes"
  ),
  make_option("--count_threshold",
    type = "integer", default = 64,
    help = "threshold for count of reads per gene to be included in plot"
  ),
  make_option("--do_pos_sp_nt_freq",
    type = "logical", default = TRUE,
    help = "do calculate the position-specific nucleotide frequency"
  ),
  make_option("--nnt_buffer",
    type = "integer", default = 25,
    help = "n nucleotides of UTR buffer to include in metagene plots"
  ),
  make_option("--nnt_gene",
    type = "integer", default = 50,
    help = "nnucleotides of gene to include in metagene plots"
  )
)

# Read in commandline arguments
# Flic: length(opt) = length(option_list) + 1 as OptionParse(add_help_option = TRUE) default is left as-is
  # this means "a standard help option should be automatically added to the OptionParser instance"
  # didn't remove that default, so length(opt) still =18
opt <- parse_args(OptionParser(option_list=option_list))
attach(opt)

print("generate_stats_figs.R running with parameters:")
opt

# Flic: checking for any NULLs passed in from other scripts wouldn't work due to NULL behaviour as this test shows:
  # length(opt) # 18, including the extra help
  # opt$features_file <- NULL
  # length(opt) # 17, setting it to NULL has dropped opt$features
# Instead need to check for any missing options.

# How best to check ^? Simple but not clever test: 
  # if(length(opt)==18){print("no dropped options")} 
# This isn't ideal for example if we then edit out the extra help test; 
# we'd be looking for 17 parameters & test would show FALSE despite script working as intended
# check by name
  # opt$features_file 
  # opt$t_rna
  # opt$codon_pos

# THOUGHT: 
# If the arguments were passed into script as NULL, they'd be immediately dropped
# Therefore if they were missing, they'd be created and set to the default?
# Therefore if they were brought in as NULL, they'd be set to NA
# So they should end up here anyway and work as expected as NA values?

### Prepare files
fid <- H5Fopen(hdFile) # Filehandle for the h5 file

# # Read in the positions of all exons/genes in GFF format and subset CDS locations
genes <- h5ls(fid,recursive = 1)$name

# Read in coding sequences
cod_seq <- readDNAStringSet(orf_fasta)

# Range of read lengths 
read_range <- MinReadLen:MaxReadLen

# Read in the positions of all exons/genes in GFF format and subset CDS locations
gff <- readGFFAsGRanges(orf_gff_file)

#####################################################################################
#####################################################################################
### Check for 3nt periodicity
print("Starting: Check for 3nt periodicity")

get_gene_datamat <- function(gene,dataset,fid) {
  # Get data matrix of read counts for gene and dataset from hd5 file fid
  data_mat <- H5Dread(H5Dopen(fid,paste0("/",gene,"/",dataset,"/reads/data"))) 
  return(data_mat)
}

get_gene_datamat_5start <- function(gene,dataset,fid,gfff,n_buffer=25,n_gene=50) {
  # get data matrix of read counts from n_buffer before start codon to n_gene after
  # for gene and dataset from hd5 file fid, using UTR5 annotations in gfff
  # if n_buffer bigger than length n_utr5, pad with zeros.
  data_mat_all <- get_gene_datamat(gene,dataset,fid)  
  n_utr5 <- width(gfff[gfff$type=="UTR5" & gfff$Name==gene])
  if(n_utr5 >= n_buffer) {
    # length n_utr5 bigger than n_buffer
    n_left5  <- n_utr5 - n_buffer + 1 # column to start from (5'end)
    zeropad5_mat <- matrix(0, 
                           nrow=nrow(data_mat_all),
                           ncol=0)
  } else {
    # length n_utr5 less than n_buffer
    n_left5  <- 1 # column to start from (5'end)
    zeropad5_mat <- matrix(0, 
                           nrow=nrow(data_mat_all),
                           ncol=n_buffer - n_utr5 )
  }
  n_right3 <- n_utr5 + n_gene # column to end with (3'end)
  data_mat_5start <- data_mat_all[,n_left5:n_right3] 
  return(  cbind(zeropad5_mat, data_mat_5start) )
}

get_gene_datamat_3end <- function(gene,dataset,fid,gfff,n_buffer=25,n_gene=50) {
  # get data matrix of read counts from n_gene before stop codon to n_buffer after
  # for gene and dataset from hd5 file fid, using UTR3 annotations in gfff
  # if n_buffer bigger than length n_utr3, pad with zeros.
  # CHECK startpos/off-by-one 
  data_mat_all <- get_gene_datamat(gene,dataset,fid)
  n_all  <- ncol(data_mat_all)
  n_utr3 <- width(gfff[gfff$type=="UTR3" & gfff$Name==gene])
  n_left5  <- n_all - n_utr3 - n_gene + 1 # column to start from (5'end)
  if(n_utr3 >= n_buffer) {
    # length n_utr3 bigger than n_buffer
    n_right3  <- n_all - n_utr3 + n_buffer # column to end with (3'end)
    zeropad3_mat <- matrix(0, 
                           nrow=nrow(data_mat_all),
                           ncol=0)
  } else {
    # length n_utr3 less than n_buffer
    n_right3 <- n_all # column to end with (3'end)
    zeropad3_mat <- matrix(0, 
                           nrow=nrow(data_mat_all),
                           ncol=n_buffer - n_utr3 )
  }
  data_mat_3end <- data_mat_all[,n_left5:n_right3] 
  return(  cbind(data_mat_3end, zeropad3_mat) )
}

tidy_datamat <- function(data_mat,startpos=1,startlen=1) {
  # CHECK startpos/off-by-one 
  positions <- startpos:(startpos+ncol(data_mat)-1)
  readlengths <- startlen:(startlen+nrow(data_mat)-1)
  data_mat %>%
    set_colnames(positions) %>%
    as_tibble %>%
    mutate(ReadLen=readlengths) %>%
    gather(-ReadLen,key = "Pos", value="Counts",convert=FALSE) %>%
    mutate(Pos=as.integer(Pos),Counts=as.integer(Counts))
}

plot_ribogrid <- function(tidymat) {
  ggplot(data=tidymat,aes(x=Pos,y=ReadLen,fill=Counts)) +
    geom_tile() +
    scale_fill_gradient("count",low = "white",high="darkblue") +
    theme(panel.grid=element_blank()) +
    labs(x="position of read 5' end", y="read length")
}

barplot_ribogrid <- function(tidymat,small_read_range=26:32) {
  ggplot(data=filter(tidymat,ReadLen %in% small_read_range),
  aes(x=Pos,y=Counts)) +
    geom_col() +
    facet_grid(ReadLen~.,scales="free_y") +
    scale_y_continuous(expand=c(0,0)) + 
    labs(x="position of read 5' end", y="count")
}

# get_gene_datamat_5start(gene="YAL003W",dataset="vignette",fid,gfff=gff) %>%
# tidy_datamat(startpos=-25,startlen=MinReadLen) %>%
#   plot_ribogrid
# get_gene_datamat_3end(gene="YAL003W",dataset="vignette",fid,gfff=gff) 


get_nt_period <- function(fid,gene,dataset,left,right){
  # previous version of script; not currently used
  data_mat <- get_gene_datamat(gene,dataset,fid) 
  pos_sum <- colSums(data_mat)  # Position-specific sum of reads of all lengths
  pos_sum <- pos_sum[left:(ncol(data_mat)-right)] # Ignore reads in parts of the buffer region defined by left/right
  return(pos_sum)
}

# NOTE: Do not use mclapply when accessing H5 data

# Get gene and position specific total counts for all read lengths
gene_poslen_counts_5start <- 
  lapply(genes,
         get_gene_datamat_5start,
         dataset=dataset,
         fid=fid,
         gfff=gff,
         n_buffer=nnt_buffer,
         n_gene=nnt_gene) %>%
  Reduce("+", .) # Reduce is a BiocGenerics, there is a purrr alternative

# ribogrid_5start 
gene_poslen_counts_5start %>%
  tidy_datamat(startpos=-nnt_buffer+1,startlen=MinReadLen) %>%
  plot_ribogrid %>%
  ggsave(filename = paste0(out_prefix,"_startcodon_ribogrid.pdf"),
      width=6,height=3)

gene_poslen_counts_5start %>%
  tidy_datamat(startpos=-nnt_buffer+1,startlen=MinReadLen) %>%
  barplot_ribogrid %>%
  ggsave(filename = paste0(out_prefix,"_startcodon_ribogridbar.pdf"),
      width=6,height=5)

gene_poslen_counts_3end <- 
  lapply(genes,
         get_gene_datamat_3end,
         dataset=dataset,
         fid=fid,
         gfff=gff,
         n_buffer=nnt_buffer,
         n_gene=nnt_gene) %>%
  Reduce("+", .) # Reduce is a BiocGenerics, there is a purrr alternative

# summarize by adding different read lengths
gene_pos_counts_5start <- gene_poslen_counts_5start %>%
  tidy_datamat(startpos=-nnt_buffer+1,startlen=MinReadLen) %>%
  group_by(Pos) %>%
  summarize(Counts=sum(Counts))

gene_pos_counts_3end <- gene_poslen_counts_3end %>%
  tidy_datamat(startpos=-nnt_gene+1,startlen=MinReadLen) %>%
  group_by(Pos) %>%
  summarize(Counts=sum(Counts))

gene_pos_counts_bothends <- bind_rows(
  gene_pos_counts_5start %>% mutate(End="5'"),
  gene_pos_counts_3end %>% mutate(End="3'") 
  ) %>%
  mutate(End = factor(End, levels = c("5'", "3'")) )

# Plot
nt_period_plot <- ggplot(gene_pos_counts_bothends, aes(x=Pos, y=Counts)) + 
  geom_line() + 
  facet_wrap(~End, scales="free") + 
  labs(x = "Nucleotide Position", y = "Read counts")

# Save plot and file
ggsave(nt_period_plot, filename = paste0(out_prefix,"_3nt_periodicity.pdf"))
write.table(gene_pos_counts_bothends,file=paste0(out_prefix,"_3nt_periodicity.tsv"),sep="\t",row=F,col=T,quote=F)

print("Completed: Check for 3nt periodicity")

#####################################################################################
#####################################################################################
### Distribution of lengths of all mapped reads
print("Starting: Distribution of lengths of all mapped reads")

# Read length-specific read counts stored as attributes of 'reads' in H5 file
gene_sp_read_len <- lapply(genes, function(x){H5Aread(H5Aopen(H5Gopen(fid,paste0("/",x,"/",dataset,"/reads")),"reads_by_len"))})

# Sum reads of each length across all genes
Counts <- colSums(matrix(unlist(gene_sp_read_len),ncol=length(read_range),byrow=T))

# Create a dataframe to store the output for plots/analyses
out_data <- data.frame(
  Length = read_range,
  Counts = Counts
)

# Plot
read_len_plot <- ggplot(out_data, aes(x=Length,y=Counts)) + 
  geom_bar(stat="identity")

# Save plot and file
ggsave(read_len_plot, filename = paste0(out_prefix,"_read_lengths.pdf"))
write.table(out_data,file=paste0(out_prefix,"_read_lengths.tsv"),sep="\t",row=F,col=T,quote=F)

print("Completed: Distribution of lengths of all mapped reads")

#####################################################################################
#####################################################################################
### Biases in nucleotide composition along mapped read lengths

print("Starting: Biases in nucleotide composition along mapped read lengths")

get_nt_read_pos <- function(fid,gene,dataset,lid,MinReadLen){
  reads_pos_len <- H5Dread(H5Dopen(fid,paste0("/",gene,"/",dataset,"/reads/data")))[lid,] # Get reads of a particular length
  reads_pos_len <- reads_pos_len[1:(length(reads_pos_len)-(lid+MinReadLen-1))]  # Ignore reads whose 5' ends map close to the end of the 3' buffer
  pos <- rep(1:length(reads_pos_len),reads_pos_len) # nt positions weighted by number of reads mapping to it
  pos_IR <- IRanges(start=pos,width=(lid+MinReadLen-1)) # Create an IRanges object for position-specific reads of a particular length
  return(pos_IR)
}

cons_mat <- function(gene,pos_IR,type="count",cframe=0,lid){
  pos_IR_frame <- pos_IR[start(pos_IR)%%3==cframe]  # Get position-specific reads of a particular length and ORF frame
  if(length(pos_IR_frame)){   
    pos_nt <- consensusMatrix(extractAt(cod_seq[[gene]],pos_IR_frame))[1:4,] # Get position-specific nucleotide counts
    if(type=="freq"){   
      pos_nt <- pos_nt/colSums(pos_nt)  # Select frequencies instead of counts
    }
  }else{
    pos_nt <- matrix(0,ncol=read_range[lid],nrow=4,dimnames = list(c("A","C","G","T"), NULL))
  }
  return(pos_nt)
}

comb_freq <- function(allfr){   
  alph <- c("A","C","G","T")
  nt_sp_freq <- c()
  for(i in alph){   
    nt_sp_freq <- rbind(nt_sp_freq,colSums(allfr[rownames(allfr)==i,])) # Get position-specific counts/freq of a nt across all genes for reads of a poarticular length and frame
  }
  nt_sp_freq <- t(nt_sp_freq/colSums(nt_sp_freq)) # Convert total counts to freq OR freq to normalized frequencies
  colnames(nt_sp_freq) <- alph
  return(nt_sp_freq)
}

if(do_pos_sp_nt_freq) {
    # This is in a conditional loop because it fails for some inputs
    # and has not been debugged.
    all_out <- c()
    for(lid in 1:length(read_range)){
        out <- lapply(genes,function(x){
            get_nt_read_pos(fid=fid,gene=as.character(x),dataset=dataset,lid=lid,MinReadLen = MinReadLen) # For each read length convert reads to IRanges
        })
        names(out) <- genes
        
        # Get position-specific nucleotide counts for reads in each frame
        fr0 <- mclapply(genes,function(gene){cons_mat(gene=gene,pos_IR=out[[gene]],cframe=0,lid=lid)},mc.cores=Ncores)
        allfr0 <- do.call(rbind,fr0)
        fr1 <- mclapply(genes,function(gene){cons_mat(gene=gene,pos_IR=out[[gene]],cframe=1,lid=lid)},mc.cores=Ncores)
        allfr1 <- do.call(rbind,fr1)
        fr2 <- mclapply(genes,function(gene){cons_mat(gene=gene,pos_IR=out[[gene]],cframe=2,lid=lid)},mc.cores=Ncores)
        allfr2 <- do.call(rbind,fr2)
        
        # Get position-specific freq for all nts
        cnt_fr0 <- signif(comb_freq(allfr0),3)
        cnt_fr1 <- signif(comb_freq(allfr1),3)
        cnt_fr2 <- signif(comb_freq(allfr2),3)
        
        output <- data.frame(rbind(cnt_fr0,cnt_fr1,cnt_fr2))
        all_out <- rbind(all_out,output)
    }
    
    # Prepare variables for output file
    Length <- unlist(lapply(read_range,function(x){rep(x,x*3)}))
    Position <- unlist(lapply(read_range,function(x){rep(1:x,3)}))
    Frame <- unlist(lapply(read_range,function(x){rep(0:2,each=x)}))
    
    all_out <- cbind(Length,Position,Frame,all_out)
    all_out[is.na(all_out)] <- 0
    
    # Save file
    write.table(all_out,file=paste0(out_prefix,"_pos_sp_nt_freq.tsv"),sep="\t",row=F,col=T,quote=F)
    
    print("Completed nucleotide composition bias table")
}

#####################################################################################
#####################################################################################
### Position specific distribution of reads

print("Starting: Position specific distribution of reads")

# Codon-specific reads for RPF datasets
get_cod_pos_reads <- function(fid,gene,dataset,left,right,MinReadLen){	
  lid <- 28-MinReadLen+1
  reads_pos <- H5Dread(H5Dopen(fid,paste0("/",gene,"/",dataset,"/reads/data"))) # Get the matrix of read counts
  reads_pos_subset <- reads_pos[,left:(dim(reads_pos)[2]-right)]  # Subset positions such that only CDS codon-mapped reads are considered
  end_reads_pos_subset <- ncol(reads_pos_subset)  # Number of columns of the subset
  
  l28 <- roll_suml(reads_pos_subset[lid,2:end_reads_pos_subset],n=3,fill=NULL)[seq(1,length(reads_pos_subset[14,2:end_reads_pos_subset]),3)]  # Map reads of length 28 to codons
  l29 <- roll_suml(reads_pos_subset[(lid+1),2:end_reads_pos_subset],n=3,fill=NULL)[seq(1,length(reads_pos_subset[15,2:end_reads_pos_subset]),3)]  # Map reads of length 29 to codons
  l30 <- roll_suml(reads_pos_subset[(lid+2),1:end_reads_pos_subset],n=3,fill=NULL)[seq(1,length(reads_pos_subset[16,1:end_reads_pos_subset]),3)]  # Map reads of length 30 to codons

  cod_sp_counts  <- l28+l29+l30   # Sum of reads of lengths 28-30 at each codon
  cod_sp_counts <- cod_sp_counts[1:(length(cod_sp_counts)-1)]
  return(cod_sp_counts)
}

# Nt-specific coverage for mRNA datasets
get_mrna_coverage <- function(fid,gene,dataset,left,right,read_range,MinReadLen,Buffer){	
  reads_pos <- H5Dread(H5Dopen(fid,paste0("/",gene,"/",dataset,"/reads/data"))) # Get the matrix of read counts
  reads_pos_subset <- reads_pos[,left:(dim(reads_pos)[2]-right)]   # Subset positions such that only CDS mapped reads are considered
  
  nt_IR_list <- lapply(read_range,function(w){IRanges(start=rep(1:ncol(reads_pos_subset),reads_pos_subset[(w-MinReadLen+1),]),width=w)})  # Create list of IRanges for position-specific reads of all length 
  nt_IR <- unlist(as(nt_IR_list,"IRangesList")) # Combine IRanges from different read lengths
  nt_cov <- coverage(nt_IR) # Estimate nt-specific coverage of mRNA reads

  # Subset coverage to only CDS
  nt_counts <- rep.int(runValue(nt_cov), runLength(nt_cov))
  if(length(nt_counts)>=(Buffer-left)){	
    nt_counts <- nt_counts[(Buffer-left):length(nt_counts)]
  }else{
    nt_counts <- 0
  }

  cds_length <- ncol(reads_pos_subset)-(Buffer-left-1)  # Length of CDS
  nt_sp_counts <- rep(0,cds_length) 

  if(length(nt_counts)<cds_length){	
    if(length(nt_counts)>0){	
      nt_sp_counts[1:length(nt_counts)] <- nt_counts
    }
  }else{
    nt_sp_counts <- nt_counts[1:cds_length]
  }
  return(nt_sp_counts)
}

# For RPF datasets, generate codon-based position-specific reads
if(rpf){
  out5p <- matrix(NA,nrow=length(genes),ncol=500) # Empty matrix to store position-specific read counts (5')
  out3p <- matrix(NA,nrow=length(genes),ncol=500) # Empty matrix to store position-specific read counts (3')

  out <- lapply(genes,function(gene){
    get_cod_pos_reads(fid=fid,gene=gene,dataset=dataset,left=(Buffer-15),right=(Buffer+11),MinReadLen = MinReadLen)}) # Get codon-based position-specific reads for each gene
  names(out) <- genes
  
  cc <- 1
  for(gene in genes){
    tmp <- out[[gene]]
    if(sum(tmp)>=64){	# Only consider genes with at least 64 mapped reads along its CDS
      tmp <- tmp/mean(tmp)
      if(length(tmp)>500){	
        out5p[cc,] <- tmp[1:500]
        out3p[cc,] <- rev(tmp)[1:500]
      }else{
        out5p[cc,1:length(tmp)] <- tmp
        out3p[cc,1:length(tmp)] <- rev(tmp)
      }
    }
    cc <- cc+1
  }

  # Estimate position-specific mean and std error of mapped read counts
  m5p <- signif(apply(out5p,2,mean,na.rm=T),4)
  m3p <- signif(apply(out3p,2,mean,na.rm=T),4)
  s5p <- signif(apply(out5p,2,function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}),4)
  s3p <- signif(apply(out3p,2,function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}),4)

  # Normalize reads to last 50 codons of the 500-codon window.
  # This allows easy comparison between datasets
  s5p <- s5p/mean(m5p[450:500])
  s3p <- s3p/mean(m3p[450:500])
  m5p <- m5p/mean(m5p[450:500])
  m3p <- m3p/mean(m3p[450:500])
  


  # Create a dataframe to store the output for plots/analyses
  out_df <- data.frame(  Position = c(1:500,0:-499),
                         Mean = c(m5p,m3p),
                         SD = c(s5p,s3p),
                         End = factor(rep(c("5'","3'"),each=500),levels=c("5'","3'"))
  )
  
  # Plot
  pos_sp_rpf_norm_reads <- ggplot(out_df, 
                                  aes(Position,Mean,col=End)) + 
    geom_line() + 
    facet_grid(~End, scales="free") + 
    guides(col=FALSE)

  # Save plot and file
  ggsave(pos_sp_rpf_norm_reads, filename = paste0(out_prefix,"_pos_sp_rpf_norm_reads.pdf"))
  write.table(out_df,file=paste0(out_prefix,"_pos_sp_rpf_norm_reads.tsv"),sep="\t",row=F,col=T,quote=F)
}

# For mRNA datasets, generate nt-based position-specific reads
if(!rpf){
  out5p <- matrix(NA,nrow=length(genes),ncol=1500)  # Empty matrix to store position-specific read counts (5')
  out3p <- matrix(NA,nrow=length(genes),ncol=1500)  # Empty matrix to store position-specific read counts (3')
  
  out <- lapply(genes,function(gene){get_mrna_coverage(fid=fid,gene=gene,dataset=dataset,left=(Buffer-49),right=(Buffer-3),MinReadLen=MinReadLen, read_range=read_range, Buffer=Buffer)})
  names(out) <- genes
  
  cc <- 1
  for(gene in genes){
    tmp <- out[[gene]]
    if(sum(tmp)>0){	# Only consider genes with at least 1 mapped read along its CDS
      tmp <- tmp/mean(tmp)
      if(length(tmp)>1500)
      {	out5p[cc,] <- tmp[1:1500]
        out3p[cc,] <- rev(tmp)[1:1500]
      }else{
        out5p[cc,1:length(tmp)] <- tmp
        out3p[cc,1:length(tmp)] <- rev(tmp)
      }
    }
    cc <- cc+1
  }
 
  # Estimate position-specific mean and std error of mapped read counts
  m5p <- signif(apply(out5p,2,mean,na.rm=T),4)
  m3p <- signif(apply(out3p,2,mean,na.rm=T),4)
  s5p <- signif(apply(out5p,2,function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}),4)
  s3p <- signif(apply(out3p,2,function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}),4)
  
  # Normalize reads to last 150 nts of the 1500-nt window.
  # This allows easy comparison between datasets
  s5p <- s5p/mean(m5p[1350:1500])
  s3p <- s3p/mean(m3p[1350:1500])
  m5p <- m5p/mean(m5p[1350:1500])
  m3p <- m3p/mean(m3p[1350:1500])
  

  # Create a dataframe to store the output for plots/analyses
  out_df <- data.frame(Position=c(1:1500,0:-1499),
                       Mean=c(m5p,m3p),
                       SD=c(s5p,s3p),
                       End=factor(rep(c("5'","3'"),each=1500),levels=c("5'","3'")))

  # Plot
  pos_sp_mrna_norm_coverage <- ggplot(out_df, aes(Position,Mean,col=End)) + 
    geom_line() + 
    facet_grid(~End, scales="free") + 
    guides(col=FALSE)
  
  # Save plot and file
  ggsave(pos_sp_mrna_norm_coverage, filename = paste0(out_prefix,"_pos_sp_mrna_norm_coverage.pdf"))
  write.table(out_df,file=paste0(out_prefix,"_pos_sp_mrna_norm_coverage.tsv"),sep="\t",row=F,col=T,quote=F)
}

print("Completed: Position specific distribution of reads")

#####################################################################################
#####################################################################################
### Calculate TPMs of genes

print("Starting: Calculate TPMs of genes")

# Read length-specific total counts stored as attributes of 'reads_total' in H5 file
get_gene_len <- function(gene,ffid=fid) {
    start_codon_pos <- H5Aread(H5Aopen(H5Gopen(ffid,paste0("/",gene,"/",dataset,"/reads")),"start_codon_pos"))[1]
    stop_codon_pos <- H5Aread(H5Aopen(H5Gopen(ffid,paste0("/",gene,"/",dataset,"/reads")),"stop_codon_pos"))[1]
    return(stop_codon_pos - start_codon_pos)
}

get_gene_reads_total <- function(gene,ffid=fid) {
    H5Aread(H5Aopen(H5Gopen(ffid,paste0("/",gene,"/",dataset,"/reads")),"reads_total"))
}

get_gene_readdensity <- function(gene,ffid=fid,buffer=50) {
    # buffer
    get_gene_reads_total(gene,ffid) / ( get_gene_len(gene,ffid) + 50 )
}

# Calculate transcripts per million (TPM)
gene_sp_reads <- sapply(genes, get_gene_reads_total, ffid=fid)
reads_per_b <- sapply(genes, get_gene_readdensity, ffid=fid)

tpms <- data.frame(ORF=genes,
                   readcount=gene_sp_reads,
                   rpb=reads_per_b,
                   tpm=reads_per_b*1e6/sum(reads_per_b) )

write.table(tpms,file=paste0(out_prefix,"_tpms.tsv"),
            sep="\t",row=F,col=T,quote=F)


#####################################################################################
#####################################################################################
### Correlations between TPMs of genes with their sequence-based features
if ( !is.na( features_file) ) {
  print("Starting: Correlations between TPMs of genes with their sequence-based features")
  
  features <- read.table(features_file,h=T)
  # features <- features[match(genes,features$ORF),] 
  # features <- features[!is.na(features$ORF),]
  # features <- merge(features,tpms,by="ORF")
  
  
  # Prepare data for plot
  # Consider only genes with at least count_threshold mapped reads
  
  plot_data <- merge(features,tpms,by="ORF") %>%
    filter(readcount >= count_threshold, !is.na(ORF)) %>%
    select(-readcount,-rpb) %>%
    gather(Feature,Value, -ORF, -tpm)
  
  features_plot <- ggplot(plot_data, aes(x=tpm, y=Value)) +
    geom_point(alpha=0.3) +
    facet_wrap(~Feature,scales="free") +
    scale_x_log10() +
    geom_smooth(method="lm") +
    xlab("TPM (transcripts per million)")
  
  # Save plot and file
  ggsave(features_plot, filename = paste0(out_prefix,"_features.pdf"))
  
  print("Completed: Correlations between TPMs of genes with their sequence-based features")
} else {
  print("Skipped: Correlations between TPMs of genes with their sequence-based features - features_file.tsv not provided")
}

#####################################################################################
#####################################################################################
### Codon-specific ribosome densities for correlations with tRNAs
if(!is.na(t_rna) & !is.na(codon_pos)) {
  print("Starting: Codon-specific ribosome densities for correlations with tRNAs")
  
  # Only for RPF datasets
  if(rpf){ 
    # This still depends on yeast-specific arguments and should be edited.
    yeast_tRNAs <- read.table(t_rna,h=T) # Read in yeast tRNA estimates
    load(codon_pos) # Position of codons in each gene (numbering ignores first 200 codons)
    # Reads in an object named "codon_pos"
    out <- lapply(genes,function(gene){
      # From "Position specific distribution of reads" plot
      get_cod_pos_reads(fid=fid,gene=gene,dataset=dataset,left=(Buffer-15),right=(Buffer+11),MinReadLen = MinReadLen)}) # Get codon-based position-specific reads for each gene
      names(out) <- genes
    
    gene_len <- sapply(out,length)  # Calculate gene length in codons
    out <- out[gene_len>201]  # Ignore genes with <=200 sense codons
    
    trim_out <- lapply(out,function(x){x[201:(length(x)-1)]}) # Trim first 200 codons and stop codon from each gene
    read_counts_trim <- sapply(trim_out,sum)  # Calculate read counts in trimmed genes
    trim_out <- trim_out[read_counts_trim>=64]  # Ignore genes with fewer than 64 mapped reads
    
    norm_out <- lapply(trim_out,function(x){x/mean(x)}) # Normalize reads in each gene by their mean
    
    # Calculate codon-specific mean ribosome-densities at A/P/E sites of the mapped reads
    a_mn <- sapply(names(codon_pos),function(codon){mean(unlist(apply(codon_pos[[codon]],1,function(a){pos <- as.numeric(a[2]);norm_out[[a[1]]][pos]})),na.rm=T)})
    p_mn <- sapply(names(codon_pos),function(codon){mean(unlist(apply(codon_pos[[codon]],1,function(a){pos <- as.numeric(a[2])+1;norm_out[[a[1]]][pos]})),na.rm=T)})
    e_mn <- sapply(names(codon_pos),function(codon){mean(unlist(apply(codon_pos[[codon]],1,function(a){pos <- as.numeric(a[2])+2;norm_out[[a[1]]][pos]})),na.rm=T)})
    
    # Sort the values
    A <- a_mn[order(names(codon_pos))]
    P <- p_mn[order(names(codon_pos))]
    E <- e_mn[order(names(codon_pos))]
    
    out_df <- cbind(yeast_tRNAs,A,P,E)
    
    # Prepare data for plot
    plot_data <- out_df %>% 
      gather(tRNA_type,tRNA_value,3:6) %>% 
      gather(Site,Ribodens,3:5)
    
    cod_dens_tRNA_plot <- ggplot(plot_data, aes(x=tRNA_value, y=Ribodens)) +
      geom_point(alpha=0.3) +
      facet_grid(Site~tRNA_type,scales="free_x") +
      geom_smooth(method = "lm") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    
    # Save plot and file
    ggsave(cod_dens_tRNA_plot, filename = paste0(out_prefix,"_codon_ribodens.pdf"))
    write.table(out_df,file=paste0(out_prefix,"_codon_ribodens.tsv"),sep="\t",row=F,col=T,quote=F)
  }
  
  print("Completed: Codon-specific ribosome densities for correlations with tRNAs")
} else {
  print("Skipped: Codon-specific ribosome densities for correlations with tRNAs - t_rna file and/or codon_pos file not provided")
}
