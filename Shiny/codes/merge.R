setwd('/Users/xtj/Documents/shahlab/shinyProject/data/')

in1<-'RPKM/F8_RPKMs.tsv'
in2<-'RPKM/bgdb.tsv'
in3<-'data_unq.tsv'

d1<-read.delim(in1)
d2<-read.table(in2)
d2<-d2[,c(1,2,5,9,10)]
names(d2)<-c()
names(d2)<-c('Year','Author','Dataset','Type','Condition')

d3<-read.table(in3)
names(d3)=c('Year','Author','Journal','NA','Dataset','NA','NA','NA','Type','Condition','geoID')
d3<-d3[,c(1,2,5,9,10)]

d4<-merge(d1, d3, by=c("Year", "Author","Dataset"))
d4<-d4[,c(1,2,3,5297,5298,4:5296)]
write.table(d4,file='RPKM/F8_RPKMs_modified.tsv',sep='\t',quote=FALSE)


d5<-merge(d1, d2, by=c("Year", "Author","Dataset"))
d5<-d5[,c(1,2,3,5297,5298,4:5296)]
write.table(d5,file='RPKM/bgdb1.tsv',sep='\t',quote=FALSE)

