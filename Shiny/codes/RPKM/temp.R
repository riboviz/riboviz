library(ggplot2)

setwd('/Users/xtj/Documents/shahlab/shinyProject/data/RPKM')
t2<-read.delim('bgdb.tsv')

dat<-t2[t2$Type=='RPF',c(5,7)]
#Solution 1 for histogram
ggplot(dat, aes(x=YAL001C)) + 
  geom_histogram(data = subset(dat,Condition == 'CHX'), fill = "red", alpha = 0.2, binwidth = 1) + 
  geom_histogram(data = subset(dat,Condition == 'FF'), fill = "green", alpha = 0.2, binwidth = 1) +
  geom_vline(xintercept = 15)

#Solution 2 for histogram
ggplot(dat, aes(x=YAL002W,color=Condition,fill=Condition)) + 
  geom_histogram(bins = 10,position='dodge',alpha=0.5) +
  theme(legend.position='top')