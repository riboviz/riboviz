a5<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b5<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf_chx_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
c5<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

a3<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b3<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf_chx_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
c3<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_rpf_ff_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

m5<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_mrna.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
mb5<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_mrna_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)


m3<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_mrna.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
mb3<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_mrna_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

thefinaldf5=data.frame(NULL)
thefinaldf3=data.frame(NULL)
for (i in 1:nrow(a5)){
print(i)

newdf5=data.frame("Year"=rep(a5[i,]$Year, 500), "Author"=rep(a5[i,]$Author, 500),"Dataset"=rep(a5[i,]$Dataset, 500),"Position"=499:0, "X5prime"=rep(1, 500), "Reads"=as.numeric(a5[i,][1, 4:503]), "SD"=rep(0, 500), "Type"=rep(100, 500), "DataType"=rep("RPF", 500))
newdf3=data.frame("Year"=rep(a3[i,]$Year, 500), "Author"=rep(a3[i,]$Author, 500),"Dataset"=rep(a3[i,]$Dataset, 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=as.numeric(a3[i,][1, 4:503]), "SD"=rep(0, 500), "Type"=rep(100, 500), "DataType"=rep("RPF", 500))


thefinaldf5=rbind(thefinaldf5, newdf5)
thefinaldf3=rbind(thefinaldf3, newdf3)
}


newdfb5=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=499:0, "X5prime"=rep(1, 500), "Reads"=as.numeric(b5[1, 2:501]), "SD"=as.numeric(b5[2, 2:501]), "Type"=rep(1, 500), "DataType"=rep("RPF", 500))
newdfc5=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=499:0, "X5prime"=rep(1, 500), "Reads"=as.numeric(c5[1, 2:501]), "SD"=as.numeric(c5[2, 2:501]), "Type"=rep(0, 500), "DataType"=rep("RPF", 500))
newdfb3=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=as.numeric(b3[1, 2:501]), "SD"=as.numeric(b3[2, 2:501]), "Type"=rep(1, 500), "DataType"=rep("RPF", 500))
newdfc3=data.frame("Year"=rep("All", 500), "Author"=rep("All", 500),"Dataset"=rep("All", 500),"Position"=499:0, "X5prime"=rep(0, 500), "Reads"=as.numeric(c3[1, 2:501]), "SD"=as.numeric(c3[2, 2:501]), "Type"=rep(0, 500), "DataType"=rep("RPF", 500))


thefinalrpf=rbind(thefinaldf5,newdfb5, newdfc5,thefinaldf3,newdfb3, newdfc3)
thefinalmrna=rbind(thefinaldf5,newdfb5, newdfc5,thefinaldf3,newdfb3, newdfc3)

thefinal=rbind(thefinalrpf, thefinalmrna)

write.table(thefinal,file=paste("/Users/carja/Dropbox/RiboViz/ProjectFolder3/VizData/F5.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")



// F5_3p_mrna.tsv				
// F5_3p_mrna_bg.tsv
a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_mrna.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_3p_mrna_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

a1=a[which( a$Year==2009 & a$Author=="Ingolia" & a$Dataset=="mRNA-rich-1"),]
a2=a[which( a$Year==2009 & a$Author=="Ingolia" & a$Dataset=="mRNA-starved-1"),]
a3=a[which( a$Year==2012 & a$Author=="Brar" & a$Dataset=="mRNA_B_tp_traditional_tc"),]
a4=a[which( a$Year==2014 & a$Author=="Pop" & a$Dataset=="WT_mrna"),]
a5=a[which( a$Year==2016 & a$Author=="Weinberg" & a$Dataset=="unselected_total_RNA"),]

dim(a1[,4:1503]) #1 1500

newdf1=data.frame("Year"=rep(a1$Year, 1500), "Author"=rep(a1$Author, 1500),"Dataset"=rep("rpf-rich-1", 1500),"Position"=1499:0, "X5prime"=rep(0, 1500), "Reads"=as.numeric(a1[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 1500))
newdf2=data.frame("Year"=rep(a2$Year, 500), "Author"=rep(a2$Author, 1500),"Dataset"=rep("rpf-starved-1", 1500),"Position"=1499:0, "X5prime"=rep(0, 1500), "Reads"=as.numeric(a2[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 1500))
newdf3=data.frame("Year"=rep(a3$Year, 1500), "Author"=rep(a3$Author, 1500),"Dataset"=rep("rpf_B_tp_traditional_tc", 1500),"Position"=1499:0, "X5prime"=rep(0, 1500), "Reads"=as.numeric(a3[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 1500))
newdf4=data.frame("Year"=rep(a4$Year, 1500), "Author"=rep(a4$Author, 1500),"Dataset"=rep("WT_rpf", 1500),"Position"=1499:0, "X5prime"=rep(0, 1500), "Reads"=as.numeric(a4[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 500))
newdf5=data.frame("Year"=rep(a5$Year, 1500), "Author"=rep(a5$Author, 1500),"Dataset"=rep("RPF", 1500),"Position"=1499:0, "X5prime"=rep(0, 1500), "Reads"=as.numeric(a5[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 1500))

newdfb=data.frame("Year"=rep("All", 1500), "Author"=rep("All", 1500),"Dataset"=rep("All", 1500),"Position"=1499:0, "X5prime"=rep(0, 1500), "Reads"=as.numeric(b[1, 1:1500]), "SD"=rep(0, 1500), "Type"=rep(1, 1500), "DataType"=rep("mRNA", 1500))

newdftotal3=rbind(newdf1, newdf2, newdf3, newdf4, newdf5, newdfb)

## now the 5 prime mRNA
##*********************
// F5_3p_mrna.tsv				
// F5_3p_mrna_bg.tsv
a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_mrna.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)
b<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_5p_mrna_bg.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

a1=a[which( a$Year==2009 & a$Author=="Ingolia" & a$Dataset=="mRNA-rich-1"),]
a2=a[which( a$Year==2009 & a$Author=="Ingolia" & a$Dataset=="mRNA-starved-1"),]
a3=a[which( a$Year==2012 & a$Author=="Brar" & a$Dataset=="mRNA_B_tp_traditional_tc"),]
a4=a[which( a$Year==2014 & a$Author=="Pop" & a$Dataset=="WT_mrna"),]
a5=a[which( a$Year==2016 & a$Author=="Weinberg" & a$Dataset=="unselected_total_RNA"),]

dim(a1[,4:1503]) #1 1500

newdf1=data.frame("Year"=rep(a1$Year, 1500), "Author"=rep(a1$Author, 1500),"Dataset"=rep("rpf-rich-1", 1500),"Position"=1499:0, "X5prime"=rep(1, 1500), "Reads"=as.numeric(a1[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 1500))
newdf2=data.frame("Year"=rep(a2$Year, 500), "Author"=rep(a2$Author, 1500),"Dataset"=rep("rpf-starved-1", 1500),"Position"=1499:0, "X5prime"=rep(1, 1500), "Reads"=as.numeric(a2[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 1500))
newdf3=data.frame("Year"=rep(a3$Year, 1500), "Author"=rep(a3$Author, 1500),"Dataset"=rep("rpf_B_tp_traditional_tc", 1500),"Position"=1499:0, "X5prime"=rep(1, 1500), "Reads"=as.numeric(a3[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 1500))
newdf4=data.frame("Year"=rep(a4$Year, 1500), "Author"=rep(a4$Author, 1500),"Dataset"=rep("WT_rpf", 1500),"Position"=1499:0, "X5prime"=rep(1, 1500), "Reads"=as.numeric(a4[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 500))
newdf5=data.frame("Year"=rep(a5$Year, 1500), "Author"=rep(a5$Author, 1500),"Dataset"=rep("RPF", 1500),"Position"=1499:0, "X5prime"=rep(1, 1500), "Reads"=as.numeric(a5[1, 4:1503]), "SD"=rep(0, 1500), "Type"=rep(100, 1500), "DataType"=rep("mRNA", 1500))

newdfb=data.frame("Year"=rep("All", 1500), "Author"=rep("All", 1500),"Dataset"=rep("All", 1500),"Position"=1499:0, "X5prime"=rep(1, 1500), "Reads"=as.numeric(b[1, 1:1500]), "SD"=rep(0, 1500), "Type"=rep(1, 1500), "DataType"=rep("mRNA", 1500))

newdftotal4=rbind(newdf1, newdf2, newdf3, newdf4, newdf5, newdfb)

newdf=rbind(newdftotal1, newdftotal2, newdftotal3, newdftotal4)


write.table(newdf,file=paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F5_Temp.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")

