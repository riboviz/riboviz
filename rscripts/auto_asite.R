library(reshape2)

### Dataframe with 3 columns: d5, length, count
### Lists the frequency of each length/d5 combo for reads corresponding to start codon
start_codon <- ...
### Dataframe with 3 columns: d3, length, count
### Lists the frequency of each length/d5 combo for reads corresponding to start codon
stop_codon <- ...

stop_codon <- acast(stop_codon, length~d3, value.var="count")
  
allowed_d5 <- seq(14,18)
allowed_d3 <- seq(8,12)
allowed_len <- seq(27,31)

mod_to_frame <- c(0,2,1)
allowed_frames = mod_to_frame[(allowed_d5 %% 3) + 1]

# asite <- data.frame(frame0 = c(15,15,15,18,18), 
#                     frame1 = c(14,14,17,17,17),
#                     frame2 = c(NA,16,16,16,NA),
#                     row.names = allowed_len)

asite <- data.frame(frame0 = rep(0, length(allowed_len)),
                    frame1 = rep(0, length(allowed_len)),
                    frame2 = rep(0, length(allowed_len)),
                    row.names = allowed_len)

# Check ambiguous rules
# check_stop inputs: read length, options for A site offset
# check_stop output: A site offset with higher frequency corresponding 3' digestion length
check_stop = function(len, frame) {
  options = allowed_d5[which(allowed_frames == frame)]
  # only consider valid d5/d3 combos
  options = options[which((len - options - 3) %in% allowed_d3)]
  if(length(options) == 1){
    return(options[1])
  }
  if(length(options) == 0){
    return(NA)
  }
  d3_options = len - options - 3
  d3_options = stop_codon[as.character(len), as.character(d3_options)]
  best_d3 = as.numeric(names(d3_options)[which.max(d3_options)])
  return(len - best_d3 - 3)
}


# Find ambiguities

for(len in allowed_len){
  for(frame in c(0,1,2)){
    asite[as.character(len), paste0('frame',frame)] = check_stop(len, frame)
  }
}

write.table(asite, file='asite_rules.tsv', quote=FALSE, sep='\t', row.names = TRUE, col.names = TRUE)

# asite = read.csv2(file='asite_rules.tsv', header = TRUE, sep = '\t', row.names = 1)
