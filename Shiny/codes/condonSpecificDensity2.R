library(rhdf5)
library(parallel)

setwd('C:/Users/txing/Desktop/ShinyProject/')
dls<-h5ls('test.h5')

myfun <- function(genename){
  group<-paste(genename,'/2016_Weinberg_RPF/reads/data',sep='')
  coord<-dim(h5read('test.h5',group))[2] #nucleotide coordinate e.g.3980 for YAL001C
  endcoord<-coord-247 #end coordinate for the last nucleotide in the ORF
  d1<-h5read('test.h5',group,index=list(14:15,(251-15):(endcoord-15))) #matrix 2 rows (ribosome protected length=28,29),columns=endcoord-251+1
  d2<-h5read('test.h5',group,index=list(16,(251-16):(endcoord-16))) #matrix 1 row (ribosome protected length=30)
  d<-rbind(d1,d2) #3X7482, matrix to calculate the ribosome density
  num_of_codons<-(endcoord-251+1)/3 #number of codons= integer
  dim(d)<-c(nrow(d)*3,num_of_codons) #integrate each 3x3 of matrix(d) into one column
  res<-lapply(data.frame(d),sum)
  
  print(genename)
  #  print(head(res,3)) #print genename and the density for first 3 codons
  H5close()
}

#mclapply(dls$name[seq(1,length(dls$name),5)],myfun)
mclapply(dls$name[seq(1,5,5)],myfun)

genename<-'YAL001C'
group<-paste(genename,'/2016_Weinberg_RPF/reads/data',sep='')
coord<-dim(h5read('test.h5',group))[2] #nucleotide coordinate e.g.3980 for YAL001C
endcoord<-coord-247 #end coordinate for the last nucleotide in the ORF
d1<-h5read('test.h5',group,index=list(14:15,(251-15):(endcoord-15))) #matrix 2 rows (ribosome protected length=28,29),columns=endcoord-251+1
d2<-h5read('test.h5',group,index=list(16,(251-16):(endcoord-16))) #matrix 1 row (ribosome protected length=30)
d<-rbind(d1,d2) #3X7482, matrix to calculate the ribosome density
num_of_codons<-(endcoord-251+1)/3 #number of codons= integer
dim(d)<-c(nrow(d)*3,num_of_codons) #integrate each 3x3 of matrix(d) into one column
res<-lapply(data.frame(d),sum)
plot(1:num_of_codons,res,type='l')
