### Dataframe with 3 columns: d5, length, count
### Lists the frequency of each length/d5 combo for reads corresponding to start codon
start_codon <- ...
### Dataframe with 3 columns: d3, length, count
### Lists the frequency of each length/d5 combo for reads corresponding to start codon
stop_codon <- ...

asite <- data.frame(frame0 = c(15,15,15,18,18), 
                    frame1 = c(14,14,17,17,17),
                    frame2 = c(NA,16,16,16,NA),
                    row.names = seq(27,31))

# Check ambiguous rules
# check_stop inputs: read length, two options for A site offset
# check_stop output: A site offset with higher frequency corresponding 3' digestion length
check_stop = function(length, option1, option2) {
  countd3_option1 = as.numeric(stop_codon[stop_codon$length == length & stop_codon$d3 == (length - option1 - 3), ]['count']) 
  countd3_option2 = as.numeric(stop_codon[stop_codon$length == length & stop_codon$d3 == (length - option2 - 3), ]['count'])
  if(countd3_option1 >= countd3_option2) {return(option1)}
  else {return(option2)}
}

asite['29','frame0'] = check_stop(29, 15, 18)
asite['30','frame0'] = check_stop(30, 15, 18)
asite['28','frame1'] = check_stop(28, 14, 17)
asite['29','frame1'] = check_stop(29, 14, 17)

write.table(asite, file='asite_rules.tsv', quote=FALSE, sep='\t', row.names = TRUE, col.names = TRUE)

# asite = read.csv2(file='asite_rules.tsv', header = TRUE, sep = '\t', row.names = 1)
