#!/usr/bin/env Rscript

# to run, call thus: Rscript --vanilla rmarkdown/ParamsTestScript.R Goose SomeValue 3 8
  # parameter1: myName (defaults in render() call DO NOT OVERWRITE commandline arguments)
  # parameter2: defaultParam (defaults in render() call DO NOT OVERWRITE commandline arguments)
  # parameter3: min_read_length (defaults in render() call DO NOT OVERWRITE commandline arguments)
  # parameter4: max_read_length (defaults in render() call DO NOT OVERWRITE commandline arguments)
  # parameter5: NEW parameter created by render() default; not fed in from commandline argument
  # parameter6: NOT FED IN VIA COMMANDLINE CALL OR RENDER ARGUMENT; created from .Rmd default parameters

args = commandArgs(trailingOnly=TRUE)


## @knitr intialargs_params
# print args
cat(args, sep = "\n")

myName <- args[1]
defaultParam <- args[2]
min_read_length <- args[3]
max_read_length <- args[4]

paste("In ParamsTestScript.R, Min_read_length is", min_read_length)
paste("In ParamsTestScript.R, Max_read_length is", max_read_length)


## @knitr check_params
# range of read lengths between parameters set in config file
read_range <- min_read_length:max_read_length

read_range # 1 line away, still part of check_params block

read_range - 2  # 2 lines away, still part of check_params block 

## @knitr junk_block_I_wont_use
print("This is a junk block I don't want to use or ever see again")

## @knitr other_block
# some other block

x <- 10
y <- 20


## @knitr makefileblock

fileConn<-file("./ParamTestFileOutput.txt")
writeLines(args, fileConn)
close(fileConn)

## @knitr metablock

rmarkdown::render(
  "rmarkdown/ParamsTest.Rmd", 
  params=list(myName="Flic", 
              defaultParam="NotDefault", 
              min_read_length=5, 
              max_read_length=25, 
              extraParam="What"
              )
  )
# parameter1: myName (defaults in render() call DO NOT OVERWRITE commandline arguments)
# parameter2: defaultParam (defaults in render() call DO NOT OVERWRITE commandline arguments)
# parameter3: min_read_length (defaults in render() call DO NOT OVERWRITE commandline arguments)
# parameter4: max_read_length (defaults in render() call DO NOT OVERWRITE commandline arguments)
# parameter5: extraParam NEW parameter created by render() default; not fed in from commandline argument
# parameter6: NOT FED IN VIA COMMANDLINE CALL OR RENDER ARGUMENT; will be created from .Rmd default parameters
