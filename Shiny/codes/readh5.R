library(rhdf5)

path<-'/data/riboseq/'
inf='2015/Sen/mRNA_WT_DED1_rep_1_15_deg/GSM1621992.h5'
head(h5ls(inf),20)
h5readAttributes(file=paste0(path,inf),name='YAL001C/2015_Sen_mRNA_WT_DED1_rep_1_15_deg/reads')
dim(h5read(inf,'YAL001C/2016_Weinberg_RPF/reads/data'))
head(h5read(inf,'YAL001C/2016_Weinberg_RPF/reads/data'))
#read submatrix in hdf5
h5read(inf,'YAL001C/2016_Weinberg_RPF/reads/data',index=list(14:16,251:260))
res<-lapply(data.frame(data),sum)
plot(1:dim(data)[2],res,type='l')


inf='/Users/xtj/Documents/shahlab/shinyProject/data/riboseq/2015/Sen/ribo_ded1-cs_rep_1_15_deg/GSM1621990.h5'
h5readAttributes(file=inf,name='YAL001C/2015_Sen_ribo_ded1-cs_rep_1_15_deg/reads')
h5read(inf,'YAL001C/2015_Sen_ribo_ded1-cs_rep_1_15_deg/reads/data')
h5read(inf,'YAL001C/2015_Sen_ribo_ded1-cs_rep_1_15_deg/reads/data',index=list(14,1:3980))
