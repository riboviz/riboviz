library(rhdf5)
library(ggplot2)
library(reshape2)

setwd('/data/riboseq/mouse/Liver/Ribosome_protected_fragment_data/Sample_1C-RP/')
inf='Sample_1C-RP.h5'

coord<-dim(h5read(inf,'ENSMUST00000000049.5/Liver_Ribosome_protected_fragment_data_Sample_1C-RP/reads/data'))[2]
d<-h5read(inf,'ENSMUST00000000049.5/Liver_Ribosome_protected_fragment_data_Sample_1C-RP/reads/data')
sumd<-lapply(data.frame(d),sum)


x<-1:coord
y<-data.frame(matrix(c(as.integer(sumd),as.integer(sumd)),ncol=2))
colnames(y)<-c('sample1','sample2')
datf <- data.frame(x, y)
dat <- melt(datf, id = 'x')


#------------------------------------------------
x <- seq(0, 4 * pi, 0.1)
n <- length(x)
y1 <- 0.5 * runif(n) + sin(x)
y2 <- 0.5 * runif(n) + cos(x) - sin(x)

df<-data.frame(x,y1,y2)
df.melted <- melt(df, id = "x")
ggplot(data = df.melted, aes(x = x, y = value, color = variable)) +
  geom_line()