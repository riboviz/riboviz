get_nt_period <- function(fid,gene,dataset,left,right)
{
	tmp <- H5Dread(H5Dopen(fid,paste("/",gene,"/",dataset,"/reads/data",sep="")))
	s <- apply(tmp,2,sum)
	s <- s[left:(dim(tmp)[2]-right)]
	return(s)
}

get_nt_read_pos <- function(fid,gene,dataset,w,left,right)
{   id <- w-14
	tmp <- H5Dread(H5Dopen(fid,paste("/",gene,"/",dataset,"/reads/data",sep="")))[id,]
	tmp <- tmp[left:(length(tmp)-right)]
	pos <- which(tmp>0)
	at <- IRanges(start=pos,width=w)

	return(at)
}   

cons_mat <- function(x,type="count",cframe=0)
{   tmp <- out[[x]][start(out[[x]])%%3==cframe]
	if(length(tmp))
	{   a <- consensusMatrix(extractAt(seq[[x]],tmp))[1:4,]
		if(type=="freq")
		{   a <- a/colSums(a)
		}
	}else{
		a <- nullmat
	}
	return(a)
}

comb_freq <- function(allfr)
{   alph <- c("A","C","G","T")
	tmp <- c()
	for(i in alph)
	{   tmp <- rbind(tmp,colSums(allfr[rownames(allfr)==i,]))
	}
	tmp <- t(tmp/colSums(tmp))
	colnames(tmp) <- alph
	return(tmp)
}

get_cod_pos <- function(fid,gene,dataset,left=left,right=right)
{	tmp <- H5Dread(H5Dopen(fid,paste("/",gene,"/",dataset,"/reads/data",sep="")))
  s <- tmp[,left:(dim(tmp)[2]-right)]
  end_s <- dim(s)[2]
	l28 <- roll_suml(s[14,2:end_s],n=3,fill=NULL)[seq(1,length(s[14,2:end_s]),3)]
	l29 <- roll_suml(s[15,2:end_s],n=3,fill=NULL)[seq(1,length(s[15,2:end_s]),3)]
	l30 <- roll_suml(s[16,1:end_s],n=3,fill=NULL)[seq(1,length(s[16,1:end_s]),3)]
	
	q <- l28+l29+l30
	q <- q[1:(length(q)-1)]
	return(q)
}

get_mrna_coverage <- function(fid,gene,dataset,left=left,right=right)
{	tmp <- H5Dread(H5Dopen(fid,paste("/",gene,"/",dataset,"/reads/data",sep="")))
	s <- tmp[,left:(dim(tmp)[2]-right)]
	a <- lapply(15:50,function(w){IRanges(start=which(s[w-14,]>0),width=w)})
	b <- RangesList(a)
	x <- IRanges(start=unlist(start(b)),width=unlist(width(b)))
	y <- coverage(x)
	z <- rep.int(runValue(y), runLength(y))
	if(length(z)>50)
	{	z <- z[51:length(z)]
	}else{
		z <- 0
	}

	l <- dim(s)[2]-50
	q <- rep(0,l)
	
	if(length(z)<l)
	{	if(length(z)>0)
		{	q[1:length(z)] <- z
		}
	}else{
		q <- z[1:l]
	}
	return(q)
}

