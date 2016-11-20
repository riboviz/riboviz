library(rhdf5)

setwd('/Users/xtj/Documents/shahlab/shinyProject/')
df<-read.table('data/data_unq.tsv')
names(df)=c('year','author','journal','NA','database','NA','NA','NA','dbtype','condition','geoID')

temp<-c("2015/Sen/ribo_ded1-cs_rep_1_15_deg", "2015/Sen/ribo_WT_DED1_rep_2_37_deg", "2015/Sen/ribo_ded1-cs_rep_2_15_deg")
input_txt<-'YAL001C'
geoid<-vector()
inf<-vector()
transtemp<-vector()
attrdat<-list()

for (i in 1:length(temp)){
  geoid[i]<-as.character(df[(df$year==unlist(strsplit(temp[i],'/'))[1])&(df$author==unlist(strsplit(temp[i],'/'))[2])&(df$database==unlist(strsplit(temp[i],'/'))[3]),11])
  inf[i]<-paste('/Users/xtj/Documents/shahlab/shinyProject/data/riboseq/',temp[i],'/',geoid[i],'.h5',sep='')
  transtemp[i]<-paste(unlist(strsplit(temp[i],'/'))[1],unlist(strsplit(temp[i],'/'))[2],unlist(strsplit(temp[i],'/'))[3],sep='_')
  attrdat[[i]]<-h5readAttributes(file=inf[i],name=paste(input_txt,'/',transtemp[i],'/reads',sep=''))$reads_by_len
}

attr_length<-(h5readAttributes(file=inf[1],name=paste(input_txt,'/',transtemp[1],'/reads',sep='')))$lengths
attr_reads_by_len<-matrix(unlist(attrdat),ncol=length(temp))
matplot(attr_length,attr_reads_by_len, type = 'l',pch=2,col =1:3)
legend("topleft", legend = temp, col=1:3,pch=1) 