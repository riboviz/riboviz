library(Rsamtools, quietly=T, warn.conflicts=F, verbose=F)
library(rtracklayer, quietly=T, warn.conflicts=F, verbose=F)
library(rhdf5, quietly=T, warn.conflicts=F, verbose=F)
library(parallel, quietly=T, warn.conflicts=F, verbose=F)
library(data.table, quietly=T, warn.conflicts=F, verbose=F)
library(RcppRoll, quietly=T, warn.conflicts=F, verbose=F)

rm(list=ls(all=T))

args=(commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
    quit(status=1)
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }   
}

source("./suppl_functions.R")

print(dataset)
print(hdfile)
print(rpf)

rpf <- as.numeric(rpf)

fid <- H5Fopen(hdfile)
genes <- h5ls(fid,recursive = 1)$name

#####################################################################################
#####################################################################################
# Figure 1 - Check 3nt periodicity of reads

out <- lapply(genes,function(gene){get_nt_period(fid=fid,gene=gene,dataset=dataset,left=226,right=225)})
names(out) <- genes

pos <- -20:25
out_df_5p <- data.frame(Position=pos)
out_df_3p <- data.frame(Position=pos)

l <- length(pos)

outmat <- lapply(out,function(x){x[1:l]})
a <- matrix(unlist(outmat),ncol=l,byrow=T)
out_df_5p$Counts <- apply(a,2,sum,na.rm=T)
out_df_5p$End <- 5

outmat <- lapply(out,function(x){rev(x)[1:l]})
a <- matrix(unlist(outmat),ncol=l,byrow=T)
out_df_3p$Counts <- rev(apply(a,2,sum,na.rm=T))
out_df_3p$End <- 3

out_df <- rbind(out_df_5p,out_df_3p)

write.table(out_df,file=paste0("F1_",dataset,".tsv"),sep="\t",row=F,col=T,quote=F)


#####################################################################################
#####################################################################################
# Figure 2 - Distribution of lengths of all mapped reads

out_df <- data.frame(Length=15:50)

out <- lapply(genes, function(x){H5Aread(H5Aopen(H5Gopen(fid,paste("/",x,"/",dataset,"/reads",sep="")),"reads_by_len"))})
a <- matrix(unlist(out),ncol=36,byrow=T)
out_df$Counts <- apply(a,2,sum,na.rm=T)

write.table(out_df,file=paste("F2_",dataset,".tsv",sep=""),sep="\t",row=F,col=T,quote=F)

#####################################################################################
#####################################################################################
# Figure 3 - Biases in nucleotide composition along mapped read lengths

if(rpf)
{ chx_bg <- read.table("./F3_RPF_CHX_BG.tsv",h=T)
  ff_bg <- read.table("./F3_RPF_FF_BG.tsv",h=T)
  bg <- rbind(chx_bg,ff_bg)

  outfile <- paste("F3_",dataset,".tsv",sep="")

  Length <- unlist(lapply(15:50,function(x){rep(x,x*3)}))
  Position <- unlist(lapply(15:50,function(x){rep(1:x,3)}))
  Frame <- unlist(lapply(15:50,function(x){rep(0:2,each=x)}))
  ASD <- CSD <- GSD <- TSD <- rep(0,length(Length))
  CHX <- 100
  blanks <- data.frame(ASD,CSD,GSD,TSD,CHX)

  all_out <- c()
  for(w in 15:50)
  {   out <- lapply(genes,function(gene){get_nt_read_pos(fid=fid,gene=gene,dataset=dataset,w=w,left=226,right=225)})

      fr0 <- lapply(1:length(genes),function(x){cons_mat(x,cframe=0)})
      allfr0 <- do.call(rbind,fr0)
      fr1 <- lapply(1:length(genes),function(x){cons_mat(x,cframe=1)})
      allfr1 <- do.call(rbind,fr1)
      fr2 <- lapply(1:length(genes),function(x){cons_mat(x,cframe=2)})
      allfr2 <- do.call(rbind,fr2)

      cnt_fr0 <- signif(comb_freq(allfr0),3)
      cnt_fr1 <- signif(comb_freq(allfr1),3)
      cnt_fr2 <- signif(comb_freq(allfr2),3)

      output <- data.frame(rbind(cnt_fr0,cnt_fr1,cnt_fr2))
	    all_out <- rbind(all_out,output)
  }
  all_out <- cbind(Length,Position,Frame,all_out,blanks)
  all_out <- rbind(all_out,bg)
  all_out[is.na(all_out)] <- 0
  write.table(all_out,file=outfile,sep="\t",row=F,col=T,quote=F)
}

#####################################################################################
#####################################################################################
# Figure 4 - Position specific distribution of reads

if(rpf)
{ bg <- read.table("./F4_RPF_BG.tsv",h=T)
  DataType <- "RPF"

out5p <- matrix(NA,nrow=length(genes),ncol=500)
out3p <- matrix(NA,nrow=length(genes),ncol=500)

outfile <- paste("F4_",dataset,".tsv",sep="")
out <- lapply(genes,function(gene){get_cod_pos(fid=fid,gene=gene,dataset=dataset,left=235,right=261)})
names(out) <- genes
cc <- 1
for(gene in genes)
{
	tmp <- out[[gene]]
	if(sum(tmp)>=64)
	{	tmp <- tmp/mean(tmp)
		if(length(tmp)>500)
		{	out5p[cc,] <- tmp[1:500]
			out3p[cc,] <- rev(tmp)[1:500]
		}else{
			out5p[cc,1:length(tmp)] <- tmp
			out3p[cc,1:length(tmp)] <- rev(tmp)
		}
	}
	cc <- cc+1
}

m5p <- signif(apply(out5p,2,mean,na.rm=T),4)
m3p <- signif(apply(out3p,2,mean,na.rm=T),4)
s5p <- signif(apply(out5p,2,function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}),4)
s3p <- signif(apply(out3p,2,function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}),4)

s5p <- s5p/mean(m5p[450:500])
s3p <- s3p/mean(m3p[450:500])
m5p <- m5p/mean(m5p[450:500])
m3p <- m3p/mean(m3p[450:500])

Position <- c(1:500,-499:0)
Mean <- c(m5p,m3p)
SD <- c(s5p,s3p)
End <- c(rep(5,500),rep(3,500))
BG <- "Data"

out_df <- data.frame(Position,Mean,SD,End,BG,DataType)
out_df <- rbind(out_df,bg)
write.table(out_df,file=outfile,sep="\t",row=F,col=T,quote=F)
}

# Figure 4 mRNA

if(!rpf)
{
  bg <- read.table("./F4_mRNA_BG.tsv",h=T)
  DataType <- "mRNA"
  
out5p <- matrix(NA,nrow=length(genes),ncol=1500)
out3p <- matrix(NA,nrow=length(genes),ncol=1500)

outfile <- paste("F4_",dataset,".tsv",sep="")
out <- lapply(genes,function(gene){get_mrna_coverage(fid=fid,gene=gene,dataset=dataset,left=201,right=247)})
names(out) <- genes
cc <- 1
for(gene in genes)
{
	tmp <- out[[gene]]
	if(sum(tmp)>0)
	{	tmp <- tmp/mean(tmp)
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

m5p <- signif(apply(out5p,2,mean,na.rm=T),4)
m3p <- signif(apply(out3p,2,mean,na.rm=T),4)
s5p <- signif(apply(out5p,2,function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}),4)
s3p <- signif(apply(out3p,2,function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}),4)

s5p <- s5p/mean(m5p[1350:1500])
s3p <- s3p/mean(m3p[1350:1500])
m5p <- m5p/mean(m5p[1350:1500])
m3p <- m3p/mean(m3p[1350:1500])

Position <- c(1:1500,-1499:0)
Mean <- c(m5p,m3p)
SD <- c(s5p,s3p)
End <- c(rep(5,1500),rep(3,1500))
BG <- "Data"

out_df <- data.frame(Position,Mean,SD,End,BG,DataType)
out_df <- rbind(out_df,bg)
write.table(out_df,file=outfile,sep="\t",row=F,col=T,quote=F)
}

#####################################################################################
#####################################################################################
# Figure 5 - Excess codon-specific ribosome densities in 5' ends of genes

if(rpf)
{ bg <- read.table("./F5_BG.tsv",h=T)
  load("./codon_pos_l250.RData")
outfile <- paste("F5_",dataset,".tsv",sep="")
out <- lapply(genes,function(gene){get_cod_pos(fid=fid,gene=gene,dataset=dataset,left=235,right=261)})
names(out) <- genes

l <- unlist(lapply(out,length))
out <- out[l>251]
s <- unlist(lapply(out,function(x){sum(x[1:(length(x)-1)])}))
out <- out[s>=64]
norm_out <- lapply(out,function(x){x/mean(x)})

a_mn_5p <- unlist(lapply(names(codon_pos),function(codon){mean(unlist(apply(codon_pos[[codon]],1,function(a){pos <- as.numeric(a[2]); pos <- pos[pos<201]; norm_out[[a[1]]][pos]})),na.rm=T)}))
a_mn_3p <- unlist(lapply(names(codon_pos),function(codon){mean(unlist(apply(codon_pos[[codon]],1,function(a){pos <- as.numeric(a[2]); pos <- pos[pos>200]; norm_out[[a[1]]][pos]})),na.rm=T)}))

a_mn_5p <- a_mn_5p[order(names(codon_pos))]+1e-3
a_mn_3p <- a_mn_3p[order(names(codon_pos))]+1e-3

Codon <- names(codon_pos)[order(names(codon_pos))]
Mean <- signif(log2(a_mn_5p/a_mn_3p),4)
SD <- 0
BG <- "Data"
out_df <- data.frame(Codon,Mean,SD,BG)
out_df <- rbind(out_df,bg)

write.table(out_df,file=outfile,sep="\t",row=F,col=T,quote=F)
}

#####################################################################################
#####################################################################################
# Figure 6 - Codon-specific ribosome densities for correlations with tRNAs

if(rpf)
{ bg <- read.table("./F6_BG.tsv",h=T)
  load("./codon_pos_i200.RData")
outfile <- paste("F6_",dataset,".tsv",sep="")
out <- lapply(genes,function(gene){get_cod_pos(fid=fid,gene=gene,dataset=dataset,left=235,right=261)})
names(out) <- genes

l <- unlist(lapply(out,length))
out <- out[l>201]
tout <- lapply(out,function(x){x[201:(length(x)-1)]})
s <- unlist(lapply(tout,function(x){sum(x[1:(length(x)-1)])}))
tout <- tout[s>=64]
norm_out <- lapply(tout,function(x){x/mean(x)})

a_mn <- unlist(lapply(names(codon_pos),function(codon){mean(unlist(apply(codon_pos[[codon]],1,function(a){pos <- as.numeric(a[2]);norm_out[[a[1]]][pos]})),na.rm=T)}))
p_mn <- unlist(lapply(names(codon_pos),function(codon){mean(unlist(apply(codon_pos[[codon]],1,function(a){pos <- as.numeric(a[2])+1;norm_out[[a[1]]][pos]})),na.rm=T)}))
e_mn <- unlist(lapply(names(codon_pos),function(codon){mean(unlist(apply(codon_pos[[codon]],1,function(a){pos <- as.numeric(a[2])+2;norm_out[[a[1]]][pos]})),na.rm=T)}))

A <- a_mn[order(names(codon_pos))]
P <- p_mn[order(names(codon_pos))]
E <- e_mn[order(names(codon_pos))]

out_df <- cbind(bg,A,P,E)
write.table(out_df,file=outfile,sep="\t",row=F,col=T,quote=F)
}

#####################################################################################
#####################################################################################
# Figure 7 - Correlations between RPKMs of genes with their sequence-based features

bg <- read.table("./F7_BG.tsv",h=T)
out <- lapply(genes, function(x){a <- h5read(fid,paste("/",x,"/",dataset,"/reads/data",sep="")); l<- dim(a)[2]-450; s <- sum(a[,226:(dim(a)[2]-225)]);c(s,l)})
a <- matrix(unlist(out),ncol=2,byrow=T)
b <- a[,1]*1e9/(sum(a[,1])*a[,2])

names(b) <- genes
outfile <- paste("F7_",dataset,".tsv",sep="")

bg$Data <- b
write.table(bg,file=outfile,sep="\t",row=F,col=T,quote=F)
