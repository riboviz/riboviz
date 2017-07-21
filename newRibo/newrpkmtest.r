library(plotly)

carrots <- data.frame(length = rnorm(100000, 6, 2))
cukes <- data.frame(length = rnorm(50000, 7, 2.5))

#Now, combine your two dataframes into one.  First make a new column in each.
carrots$veg <- 'carrot'
cukes$veg <- 'cuke'

#and combine into your new data frame vegLengths
vegLengths <- rbind(carrots, cukes)

#now make your lovely plot
p <- ggplot(vegLengths, aes(length, fill = veg)) + geom_density(alpha = 0.2)

p <- ggplotly(p)

####################################################################
load("/home/txing/RiboTest/newRibo/riboViz.RData")
library(ggplot2);library(reshape2);library(tidyr)

dbslash<-c("2016/Weinberg/unselected_total_RNA","2016/Weinberg/RPF","2015/Sen/ribo_ded1-cs_rep_1_15_deg")
tempstr<-sapply(dbslash,function(x){strsplit(x,'/')})
rpkm<-sapply(tempstr, function(x){d1[d1$Year==x[1] & d1$Author==x[2] & d1$Dataset==x[3],6:ncol(d1)]}) # matrix 5293 X 3
#rownames(rpkm) <- c()

df_rpkm<-data.frame(rpkm)
tmp_df <- df_rpkm %>% gather(dbname,RPKM,1:3)
tmp_df[2] <- lapply(tmp_df[2], as.numeric)   # Convert the 2nd column to numeric
p <- ggplot(tmp_df, aes(x=RPKM, fill = dbname)) + geom_density(alpha = 0.2) + scale_x_log10()

# a<-data.frame(rpkmvalue=matrix(rpkm[,1],ncol=1))
# a$dbname<-colnames(rpkm)[1]
# 
# b<-data.frame(rpkmvalue=matrix(rpkm[,2],ncol=1))
# b$dbname<-colnames(rpkm)[2]
# 
# mm<-rbind(a,b)
# p <- ggplot(mm, aes(rpkmvalue, fill = dbname)) + geom_density(alpha = 0.2)
