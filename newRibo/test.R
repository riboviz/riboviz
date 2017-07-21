# library(plotly)
# p1 <- plot_ly(economics, x = ~date, y = ~unemploy) %>%
#   add_lines(name = ~"unemploy")
# p2 <- plot_ly(economics, x = ~date, y = ~uempmed) %>%
#   add_lines(name = ~"uempmed")
# p <- subplot(p1, p2)

a->c("2016/Weinberg/unselected_total_RNA","2016/Weinberg/rRNA_depleted_using_RiboMinus","2015/Sen/mRNA_WT_TIF1_rep_1_30_deg")
tempstr <- sapply(a,function(x){strsplit(x,'/')})
cord<- 3980
geoid <- as.character(sapply(tempstr, function(x){df[(df$year==x[1])&(df$author==x[2])&(df$database==x[3]),11]}))
inf <- as.character(mapply(function(x,y){paste0(path,x,'/',y,'.h5')}, a,geoid))
transdb <- as.character(sapply(tempstr, function(x){paste(x[1],x[2],x[3],sep="_")}))
attrdat <- as.numeric(mapply(function(x,y){h5readAttributes(file=x,name=paste0("YAL001C",'/',y,'/reads'))$reads_by_len},inf,transdb))

d <- mapply(function(x,y){return(data.frame(h5read(x,paste0("YAL001C",'/',y,'/reads/data'))))},inf,transdb, SIMPLIFY = FALSE)
dd <- mapply(function(x,y){return(data.frame(h5read(x,paste0("YAL001C",'/',y,'/reads/data'))[,251:(coord-247)]))}, inf,transdb, SIMPLIFY = FALSE)
sumd <- lapply(dd, function(x){lapply(x,sum)})

cd <- lapply(d,function(x){return(data.frame(rbindlist(list(x[14:15,(251-15):(coord-247-15)],x[16,(251-16):(coord-247-16)]))))})
num_of_codons <- (coord-247-251+1)/3
newcd <- lapply(cd, function(x){lapply(data.frame(matrix(unlist(x),nrow=nrow(x)*3,ncol=num_of_codons)),sum)})


