library(rhdf5)

setwd('/data/riboseq/mouse/Liver/Ribosome_protected_fragment_data/Sample_1C-RP/')
inf='Sample_1C-RP.h5'
head(h5ls(inf),20)

h5readAttributes(file=inf,name="ENSMUST00000000049.5/Liver_Ribosome_protected_fragment_data_Sample_1C-RP/reads")
dim(h5read(inf,'ENSMUST00000000049.5/Liver_Ribosome_protected_fragment_data_Sample_1C-RP/reads/data'))
head(h5read(inf,'ENSMUST00000000049.5/Liver_Ribosome_protected_fragment_data_Sample_1C-RP/reads/data'))
d<-h5read(inf,'ENSMUST00000000049.5/Liver_Ribosome_protected_fragment_data_Sample_1C-RP/reads/data')
sumd<-lapply(data.frame(d),sum)

#-----------------------------------------
setwd('/data/riboseq/mouse/Liver/Total_RNA_data/Sample_1C-T/')
inf='Sample_1C-T.h5'
head(h5ls(inf),20)

h5readAttributes(file=inf,name="ENSMUST00000025253.11/Liver_Total_RNA_data_Sample_1C-T/reads")
dim(h5read(inf,'ENSMUST00000025253.11/Liver_Total_RNA_data_Sample_1C-T/reads/data'))
head(h5read(inf,'ENSMUST00000025253.11/Liver_Total_RNA_data_Sample_1C-T/reads/data'))
coord<-dim(h5read(inf,'ENSMUST00000025253.11/Liver_Total_RNA_data_Sample_1C-T/reads/data'))[2]
d<-h5read(inf,'ENSMUST00000025253.11/Liver_Total_RNA_data_Sample_1C-T/reads/data',index=list(1:36,251:(coord-247)))
sumd<-lapply(data.frame(d),sum)








#read submatrix in hdf5
h5read(inf,'YAL001C/2016_Weinberg_RPF/reads/data',index=list(14:16,251:260))
res<-lapply(data.frame(data),sum)
plot(1:dim(data)[2],res,type='l')


#read lncrna file
setwd('/data/riboseq/mouse/Muscle/Ribosome_protected_fragment_data/Sample_M-RP-C/') #genename=10006...
inf='Sample_M-RP-C_lncrna.h5'
head(h5ls(inf),20)
h5readAttributes(file=inf,name="ENSMUST00000000003.13/Muscle_Ribosome_protected_fragment_data_Sample_M-RP-C/reads")
dim(h5read(inf,'ENSMUST00000000003.13/Muscle_Ribosome_protected_fragment_data_Sample_M-RP-C/reads/data'))
head(h5read(inf,'ENSMUST00000000003.13/Muscle_Ribosome_protected_fragment_data_Sample_M-RP-C/reads/data'))

#----------------h5readAttributes example
# $`3p_utr_length`
# [1] 237
# 
# $`5p_utr_length`
# [1] 140
# 
# $buffer_left
# [1] 250
# 
# $buffer_right
# [1] 247
# 
# $lengths
# [1] 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
# 
# $reads_by_len
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# 
# $reads_total
# [1] 0
# 
# $start_codon_pos
# [1] 391 392 393
# 
# $stop_codon_pos
# [1] 913 914 915