library(rhdf5)

setwd('/Users/xtj/Documents/shahlab/shinyProject/data/')
dls<-h5ls('test.h5')
i=-4
#while(i<(length(dls$name)-5)){ #pick up genenames
while(i<1){
  i=i+5
  genename<-dls$name[i] #e.g.'YAL001C'
  group<-paste(genename,'/2016_Weinberg_RPF/reads/data',sep='')
  coord<-dim(h5read('test.h5',group))[2] #nucleotide coordinate e.g.3980 for YAL001C
  endcoord<-coord-247 #end coordinate for the last nucleotide in the ORF
  d1<-h5read('test.h5',group,index=list(14:15,(251-15):(endcoord-15))) #matrix 3 rows (ribosome protected length=28,29,30),columns=endcoord-251+1
  d2<-h5read('test.h5',group,index=list(16,(251-16):(endcoord-16)))
  d<-rbind(d1,d2) #3X3483, matrix to calculate the ribosome density
  print(dim(d))
  num_of_codons<-(endcoord-251+1)/3 #number of codons= integer
  ii=-2
  v<-c() #vector storing the codon-specific density for one gene
  while ((ii+2)<(endcoord-251+1)){
    ii=ii+3
    v<-c(v,sum(d[,ii:(ii+2)])) #e.g.sum(d[,1:3])=ribosome density for the first 3 nucleotides (1st codon) of protected length=28/29/30
  }
  print(genename)
  print(head(v,10)) #print genename and the density for first 10 codons
  H5close()
}


