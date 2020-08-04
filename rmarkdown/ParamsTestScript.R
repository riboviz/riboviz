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

