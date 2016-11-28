

## Process Figure 1.
##******************

prime5<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F1_5p.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
prime3<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F1_3p.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)


for (i in 1:nrow(prime5)){
print(i)
newdf5=data.frame("Year"=rep(prime5[i,]$Year, 46), "Author"=rep(prime5[i,]$Author, 46),"Dataset"=rep(prime5[i,]$Dataset, 46),"position"=(-20):25, "count"=as.integer(prime5[i,][1, 4:49]), "end"=rep(5, 46))
newdf3=data.frame("Year"=rep(prime3[i,]$Year, 46), "Author"=rep(prime3[i,]$Author, 46),"Dataset"=rep(prime3[i,]$Dataset, 46),"position"=(-20):25, "count"=as.integer(prime3[i,][1, (ncol(prime3)-45):ncol(prime3)]), "end"=rep(3, 46))
newdf=rbind(newdf5, newdf3)
write.table(newdf,file=paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/F1_Year_", prime5[i,]$Year, "_Author_", prime5[i,]$Author,"_Dastaset_", prime5[i,]$Dataset,"_All.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")
}

## Process Figure 2.
##******************

f2<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F2_dist.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

for (i in 1:nrow(f2)){
print(i)
newdf=data.frame("Year"=rep(f2[i,]$Year, 36), "Author"=rep(f2[i,]$Author, 36),"Dataset"=rep(f2[i,]$Dataset, 36),"key"=15:50, "value"=as.integer(f2[i,][1, 4:39]))

write.table(newdf,file=paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/F2_Year_", prime5[i,]$Year, "_Author_", prime5[i,]$Author,"_Dastaset_", prime5[i,]$Dataset,"_All.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")
}

## Process Figure 3.
##******************

a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F3_posfreq.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F3_posfreq_rpf_chx_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
c<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F3_posfreq_rpf_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

for (i in 103:nrow(f2)){
print(i) 
newdf=data.frame("Year"=rep(f2[i,]$Year, 36), "Author"=rep(f2[i,]$Author, 36),"Dataset"=rep(f2[i,]$Dataset, 36),"key"=15:50, "value"=as.integer(f2[i,][1, 4:39]))
year=newdf$Year[1]
author=as.character(newdf$Author[1])
dataset=as.character(newdf$Dataset[1])

a1=a[which(a$Year==year & a$Author==author & a$Dataset==dataset),]
if (nrow(a1)>0 & !is.na(a1[1,5])){
if (length(which(is.na(a1[1,])))==0){
Position=as.numeric(a[1,5:ncol(a1)])
}else {Position=as.numeric(a[1,5:(which(is.na(a1[1,]))[1]-1)])}
thevec1=rep(15, 15*3)
for (i in 16:Position[length(Position)]){
	thevec1=c(thevec1, rep(i,i*3))
}

if (length(thevec1)!=length(Position))
{
new=thevec1[length(thevec1)]-1
thevec1=rep(15, 15*3)
for (i in 16:new){
	thevec1=c(thevec1, rep(i,i*3))
}
Position=Position[1:length(thevec1)]
}

Frame=as.numeric(a[2,5:(length(thevec1)+5-1)])
A=as.numeric(a1[1,5:(length(thevec1)+5-1)])
C=as.numeric(a1[2,5:(length(thevec1)+5-1)])
G=as.numeric(a1[3,5:(length(thevec1)+5-1)])
T=as.numeric(a1[4,5:(length(thevec1)+5-1)])
Length=thevec1
CHX=rep(100, length(A))
newdf=data.frame("Length"=Length, "Position"=Position, "Frame"=Frame, "A"=A, "C"=C, "G"=G, "T"=T, "ASD"=rep(0, length(Position)), "CSD"=rep(0, length(Position)), "GSD"=rep(0, length(Position)), "TSD"=rep(0, length(Position)), "CHX"=CHX)

Ab=as.numeric(b[1,3:(length(thevec1)+5-1-2)])
AbSD=as.numeric(b[2,3:(length(thevec1)+5-1-2)])
Cb=as.numeric(b[3,3:(length(thevec1)+5-1-2)])
CbSD=as.numeric(b[4,3:(length(thevec1)+5-1-2)])
Gb=as.numeric(b[5,3:(length(thevec1)+5-1-2)])
GbSD=as.numeric(b[6,3:(length(thevec1)+5-1-2)])
Tb=as.numeric(b[7,3:(length(thevec1)+5-1-2)])
TbSD=as.numeric(b[8,3:(length(thevec1)+5-1-2)])
CHXb=rep(1, length(Ab))
newdfb=data.frame("Length"=Length, "Position"=Position, "Frame"=Frame, "A"=Ab, "C"=Cb, "G"=Gb, "T"=Tb, "ASD"=AbSD, "CSD"=CbSD, "GSD"=GbSD, "TSD"=TbSD,"CHX"=CHXb)
Ac=as.numeric(c[1,3:(length(thevec1)+5-1-2)])
AcSD=as.numeric(c[2,3:(length(thevec1)+5-1-2)])
Cc=as.numeric(c[3,3:(length(thevec1)+5-1-2)])
CcSD=as.numeric(c[4,3:(length(thevec1)+5-1-2)])
Gc=as.numeric(c[5,3:(length(thevec1)+5-1-2)])
GcSD=as.numeric(c[6,3:(length(thevec1)+5-1-2)])
Tc=as.numeric(c[7,3:(length(thevec1)+5-1-2)])
TcSD=as.numeric(c[8,3:(length(thevec1)+5-1-2)])
CHXc=rep(0, length(Ac))
newdfc=data.frame("Length"=Length, "Position"=Position, "Frame"=Frame, "A"=Ac, "C"=Cc, "G"=Gc, "T"=Tc, "ASD"=AcSD, "CSD"=CcSD, "GSD"=GcSD, "TSD"=TcSD, "CHX"=CHXc)

newdffinal=rbind(newdf, newdfb, newdfc)
write.table(newdffinal,file=paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/F3_Temp_Year_",a1$Year[1],"_Author_",a1$Author[1],"_Dataset_", a1$Dataset[1],"_data.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")
}
}


## Process Figure 5.
##******************

// F5_3p_mrna.tsv				
// F5_3p_mrna_bg.tsv			
// F5_3p_rpf.tsv				
// F5_3p_rpf_chx_bg.tsv			
// F5_3p_rpf_ff_bg.tsv			
// F5_5p_mrna.tsv		
// F5_5p_mrna_bg.tsv			
// F5_5p_rpf.tsv				
// F5_5p_rpf_chx_bg.tsv		
// F5_5p_rpf_ff_bg.tsv
// "Position"	"X5prime"	"Reads"	"SD"	"Type"	"DataType"



for (i in 1:nrow(f2)){
a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_rpf.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_rpf_chx_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
c<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_rpf_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

print(i) 
newdf=data.frame("Year"=rep(f2[i,]$Year, 36), "Author"=rep(f2[i,]$Author, 36),"Dataset"=rep(f2[i,]$Dataset, 36),"key"=15:50, "value"=as.integer(f2[i,][1, 4:39]))
year=newdf$Year[1]
author=as.character(newdf$Author[1])
dataset=as.character(newdf$Dataset[1])

a1=a[which(a$Year==year & a$Author==author & a$Dataset==dataset),]

dim(a1[,4:503]) #1 500
if (nrow(a1)>0){
newdf1=data.frame("Year"=rep(a1$Year, 500), "Author"=rep(a1$Author, 500),"Dataset"=rep(a1$Dataset, 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=as.numeric(a1[1, 4:503]), "SD"=rep(0, 500), "Type"=rep(100, 500), "DataType"=rep("RPF", 500))

newdfb=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=as.numeric(b[1, 2:501]), "SD"=as.numeric(b[2, 2:501]), "Type"=rep(1, 500), "DataType"=rep("RPF", 500))
newdfc=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=as.numeric(c[1, 2:501]), "SD"=as.numeric(c[2, 2:501]), "Type"=rep(0, 500), "DataType"=rep("RPF", 500))

newdftotal1=rbind(newdf1, newdfb, newdfc)


## now the 5 prime RPF
##********************
a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf_chx_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
c<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

a1=a[which(a$Year==year & a$Author==author & a$Dataset==dataset),]

dim(a1[,4:503]) #1 500

newdf1=data.frame("Year"=rep(a1$Year, 500), "Author"=rep(a1$Author, 500),"Dataset"=rep(a1$Dataset, 500),"Position"=0:499, "X5prime"=rep(1, 500), "Reads"=as.numeric(a1[1, 4:503]), "SD"=rep(0, 500), "Type"=rep(100, 500), "DataType"=rep("RPF", 500))


newdfb=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=0:499, "X5prime"=rep(1, 500), "Reads"=as.numeric(b[1, 2:501]), "SD"=as.numeric(b[2, 2:501]), "Type"=rep(1, 500), "DataType"=rep("RPF", 500))
newdfc=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=0:499, "X5prime"=rep(1, 500), "Reads"=as.numeric(c[1, 2:501]), "SD"=as.numeric(c[2, 2:501]), "Type"=rep(0, 500), "DataType"=rep("RPF", 500))

newdftotal2=rbind(newdf1, newdfb, newdfc)


newdf=rbind(newdftotal1, newdftotal2)

write.table(newdf,file=paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/F5_Year_",a1$Year[1],"_Author_",a1$Author[1],"_Dataset_", a1$Dataset[1],"_data.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")
}
}


##now the same for the mRNA
##*************************
// F5_3p_mrna.tsv				
// F5_3p_mrna_bg.tsv

for (i in 1:nrow(f2)){
a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_mrna.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_mrna_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
#c<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_rpf_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

print(i) 
newdf=data.frame("Year"=rep(f2[i,]$Year, 36), "Author"=rep(f2[i,]$Author, 36),"Dataset"=rep(f2[i,]$Dataset, 36),"key"=15:50, "value"=as.integer(f2[i,][1, 4:39]))
year=newdf$Year[1]
author=as.character(newdf$Author[1])
dataset=as.character(newdf$Dataset[1])

a1=a[which(a$Year==year & a$Author==author & a$Dataset==dataset),]

dim(a1[,4:503]) #1 500
if (nrow(a1)>0){
newdf1=data.frame("Year"=rep(a1$Year, 500), "Author"=rep(a1$Author, 500),"Dataset"=rep(a1$Dataset, 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=as.numeric(a1[1, 4:503]), "SD"=rep(0, 500), "Type"=rep(100, 500), "DataType"=rep("mRNA", 500))

newdfb=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=as.numeric(b[1, 2:501]), "SD"=as.numeric(b[2, 2:501]), "Type"=rep(1, 500), "DataType"=rep("mRNA", 500))
newdfc=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=rep(0, 500), "SD"=rep(0, 500), "Type"=rep(0, 500), "DataType"=rep("mRNA", 500))

newdftotal1=rbind(newdf1, newdfb, newdfc)


## now the 5 prime 
##********************
a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_mrna.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_mrna_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
#c<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

a1=a[which(a$Year==year & a$Author==author & a$Dataset==dataset),]

dim(a1[,4:503]) #1 500

newdf1=data.frame("Year"=rep(a1$Year, 500), "Author"=rep(a1$Author, 500),"Dataset"=rep(a1$Dataset, 500),"Position"=0:499, "X5prime"=rep(1, 500), "Reads"=as.numeric(a1[1, 4:503]), "SD"=rep(0, 500), "Type"=rep(100, 500), "DataType"=rep("mRNA", 500))


newdfb=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=0:499, "X5prime"=rep(1, 500), "Reads"=as.numeric(b[1, 2:501]), "SD"=as.numeric(b[2, 2:501]), "Type"=rep(1, 500), "DataType"=rep("mRNA", 500))
newdfc=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=0:499, "X5prime"=rep(1, 500), "Reads"=rep(0, 500), "SD"=rep(0, 500), "Type"=rep(0, 500), "DataType"=rep("mRNA", 500))

newdftotal2=rbind(newdf1, newdfb, newdfc)


newdf=rbind(newdftotal1, newdftotal2)

write.table(newdf,file=paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/F5_Year_",a1$Year[1],"_Author_",a1$Author[1],"_Dataset_", a1$Dataset[1],"_data.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")
}
}

## Process Figure 6.
##******************



for (i in 1:nrow(f2)){

a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F6_excess5pcodribodens.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F6_excess5pcodribodens_chx_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
c<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F6_excess5pcodribodens_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

print(i) 
newdf=data.frame("Year"=rep(f2[i,]$Year, 36), "Author"=rep(f2[i,]$Author, 36),"Dataset"=rep(f2[i,]$Dataset, 36),"key"=15:50, "value"=as.integer(f2[i,][1, 4:39]))
year=newdf$Year[1]
author=as.character(newdf$Author[1])
dataset=as.character(newdf$Dataset[1])


a1=a[which(a$Year==year & a$Author==author & a$Dataset==dataset),]
dim(a1[,4:64]) #1 61
if (nrow(a1)>0){
newdf1=data.frame("Year"=rep(a1$Year, 61), "Author"=rep(a1$Author, 61),"Dataset"=rep(a1$Dataset, 61),"Codon"=names(a1[,4:64]), "excess_reads"=as.numeric(a1[,4:64]), "CHX"=rep("Data", 61), "SD"=rep(0, 61))

newdfb=data.frame("Year"=rep("All", 61), "Author"=rep("All", 61),"Dataset"=rep("All", 61),"Codon"=names(b[,2:62]), "excess_reads"=as.numeric(b[1,2:62]), "CHX"=rep("CHX", 61), "SD"=as.numeric(b[2,2:62]))
newdfc=data.frame("Year"=rep("All", 61), "Author"=rep("All", 61),"Dataset"=rep("All", 61),"Codon"=names(c[,2:62]), "excess_reads"=as.numeric(c[1,2:62]), "CHX"=rep("FF", 61), "SD"=as.numeric(c[2,2:62]))

newdftotal=rbind(newdf1,  newdfb, newdfc)
write.table(newdftotal,file=paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/F6_Year_",a1$Year[1],"_Author_",a1$Author[1],"_Dataset_", a1$Dataset[1],"_data.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")
}
}



## Process Figure 8.
##******************

#F7_codribodens.tsv			
#F7_corr.tsv				
#F7_corr_chx_bg.tsv			
#F7_corr_ff_bg.tsv			
#F7_tRNA.tsv
#F8_corr.tsv
#F8_features.tsv
#F8_RPKMs.tsv

a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F8_RPKMs.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F8_features.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
c<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F8_corr.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

for (i in 1:nrow(f2)){
print(i) 
newdf=data.frame("Year"=rep(f2[i,]$Year, 36), "Author"=rep(f2[i,]$Author, 36),"Dataset"=rep(f2[i,]$Dataset, 36),"key"=15:50, "value"=as.integer(f2[i,][1, 4:39]))
year=newdf$Year[1]
author=as.character(newdf$Author[1])
dataset=as.character(newdf$Dataset[1])


a1=a[which(a$Year==year & a$Author==author & a$Dataset==dataset),]
dim(a1[,4:64]) #1 61
if (nrow(a1)>0){
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F8_features.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b$RPF=as.numeric(a1[4:length(a1)])
b=b[-which(b$RPF==0),]
write.table(b,file=paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/F8_Year_",a1$Year[1],"_Author_",a1$Author[1],"_Dataset_", a1$Dataset[1],"_data.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")

}
}


## Process Figure 7.
##******************



################### FIgure 7
#F7_codribodens.tsv			
#F7_corr.tsv				
#F7_corr_chx_bg.tsv			
#F7_corr_ff_bg.tsv			
#F7_tRNA.tsv

a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F7_tRNA.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
a2<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F7_corr_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

a3<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F7_corr_chx_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
a4<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F7_corr.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
a5<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F7_codribodens.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)


##need to now figure out hw to create the files I need 
a51=a5[which( a5$Year==2009 & a5$Author=="Ingolia" & a5$Dataset=="rpf-rich-1"),]
Codon=names(a51[5:length(a51)])
A=as.numeric(a51[1, 5:length(a51)])
P=as.numeric(a51[2, 5:length(a51)])
E=as.numeric(a51[3, 5:length(a51)])
AA=a[,1]
tRNA=as.numeric(a[,3])
tAI=as.numeric(a[,4])
Microarray =as.numeric(a[,5])
RNA.seq=as.numeric(a[,6])
newdf=data.frame("AA"=AA, "Codon"=Codon, "A"=A, "P"=P, "E"=E, "tRNA"=tRNA, "tAI"=tAI, "Microarray"=Microarray, "RNAseq"=RNA.seq)
write.table(newdf,file=paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F7_Temp_Year_",a51$Year[1],"_Author_",a51$Author[1],"_Dataset_", a51$Dataset[1],"_data.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")

a51=a5[which( a5$Year==2009 & a5$Author=="Ingolia" & a5$Dataset=="rpf-starved-1"),]
Codon=names(a51[5:length(a51)])
A=as.numeric(a51[1, 5:length(a51)])
P=as.numeric(a51[2, 5:length(a51)])
E=as.numeric(a51[3, 5:length(a51)])
AA=a[,1]
tRNA=as.numeric(a[,3])
tAI=as.numeric(a[,4])
Microarray =as.numeric(a[,5])
RNA.seq=as.numeric(a[,6])
newdf=data.frame("AA"=AA, "Codon"=Codon, "A"=A, "P"=P, "E"=E, "tRNA"=tRNA, "tAI"=tAI, "Microarray"=Microarray, "RNAseq"=RNA.seq)
write.table(newdf,file=paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F7_Temp_Year_",a51$Year[1],"_Author_",a51$Author[1],"_Dataset_", a51$Dataset[1],"_data.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")

a51=a5[which( a5$Year==2016 & a5$Author=="Weinberg" & a5$Dataset=="RPF"),]
Codon=names(a51[5:length(a51)])
A=as.numeric(a51[1, 5:length(a51)])
P=as.numeric(a51[2, 5:length(a51)])
E=as.numeric(a51[3, 5:length(a51)])
AA=a[,1]
tRNA=as.numeric(a[,3])
tAI=as.numeric(a[,4])
Microarray =as.numeric(a[,5])
RNA.seq=as.numeric(a[,6])
newdf=data.frame("AA"=AA, "Codon"=Codon, "A"=A, "P"=P, "E"=E, "tRNA"=tRNA, "tAI"=tAI, "Microarray"=Microarray, "RNAseq"=RNA.seq)
write.table(newdf,file=paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F7_Temp_Year_",a51$Year[1],"_Author_",a51$Author[1],"_Dataset_", a51$Dataset[1],"_data.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")


############
f=f2[,1:3]
f$mRNA=rep(0, nrow(f))
for (i in 1:nrow(f)){
print (i)
if (length(grep("RNA", f[i,3]))>0){
f$mRNA[i]=1}
}
f[342,4]=1
write.table(f,file=paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/AllData.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")



a<-read.table(paste("/Users/carja/Dropbox/RiboViz/RibovizPush/Data/F3_Temp_Year_2014_Author_Pop_Dataset_WT_rpf_data.tsv",sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)