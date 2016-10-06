a<-read.table(paste("/Users/carja/Dropbox/RiboViz/ProjectFolder2/VizData/F2_dist.tsv", sep=""), header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE)

thefinaldf=data.frame(NULL)

for (i in 1:nrow(a)){
print(i)
newdf=data.frame("Year"=rep(a[i,]$Year, 36), "Author"=rep(a[i,]$Author, 36),"Dataset"=rep(a[i,]$Dataset, 36),"key"=15:50, "value"=as.integer(a[i,][1, 4:39]))
thefinaldf=rbind(thefinaldf, newdf)
}
write.table(thefinaldf,file=paste("/Users/carja/Dropbox/RiboViz/ProjectFolder3/VizData/F2.tsv",sep=""), row.names=FALSE, append=FALSE, sep="\t")



