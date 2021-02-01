# testthat tests for bam_to_h5.R
#
# This assumes the following files are in the path:
#
# rscripts/bam_to_h5.R
# vignette/input/yeast_YAL_CDS_w_250utrs.gff3
# vignette/output/WTnone/WTnone.bam
#
# The test runs bam_to_h5.R using the GFF and BAM file then validates
# the .h5 file created based upon its expected qualities given those
# of the input GFF and BAM files and the bam_to_h5.R command-line
# parameters.
#
# Note: for rapid development this does not currently run bam_to_h5.R
# but validates vignette/output/WTnone/WTnone.h5. TODO remove this
# comment when development is complete.
#
# To run interactively:
#
# test_file("rscripts/tests/testthat/test_bam_to_h5.R")
#
# To run from console:
#
# Rscript -e 'library(testthat); test_file("rscripts/tests/testthat/test_bam_to_h5.R")'

suppressMessages(library(glue, quietly = T))
suppressMessages(library(here, quietly = T))
suppressMessages(library(testthat, quietly = T))
suppressMessages(library(withr, quietly = T))
suppressMessages(library(GenomicAlignments, quietly = T))
suppressMessages(library(Rsamtools, quietly = T))

source(here::here("rscripts", "read_count_functions.R"))
print(here())
bam_to_h5 <- here::here("rscripts/bam_to_h5.R")
print(bam_to_h5)
gff_file <- here::here("vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
print(gff_file)
bam_file <- here::here("vignette/output/WTnone/WTnone.bam") # TODO remove
bam_file <- here::here("WTnone.bam") # TODO remove
print(bam_file)
h5_file <- here::here("test.h5")
h5_file <- here::here("vignette/output/WTnone/WTnone.h5") # TODO remove
h5_file <- here::here("WTnone.h5") # TODO remove
print(h5_file)

context("test_bam_to_h5.R")

delete_file <- function(file_name) {
  print("Invoking delete_file fixture...")
  if (file.exists(file_name)) {
    # file.remove(file_name) # TODO uncomment
  }
}

test_that("Run bam_to_h5.R and validate H5 file", {
  withr::defer(delete_file(h5_file))

  expect_equal(0, 0, info = "Example assertion") # TODO remove

  min_read_length <- 10
  max_read_length <- 50
  buffer <- 250
  primary_id <- "Name"
  secondary_id <- "NULL"
  dataset <- "vignette"
  is_riboviz_gff <- TRUE
  stop_in_cds <- FALSE

  num_read_counts <- max_read_length - min_read_length + 1

  bam_to_h5_cmd_template <- "Rscript --vanilla {bam_to_h5} --num-processes=1 --min-read-length={min_read_length} --max-read-length={max_read_length} --buffer={buffer} --primary-id={primary_id} --secondary-id={secondary_id} --dataset={dataset} --bam-file={bam_file} --hd-file={h5_file} --orf-gff-file={gff_file} --is-riboviz-gff={is_riboviz_gff} --stop-in-cds={stop_in_cds}"
  print(bam_to_h5_cmd_template)
  cmd <- glue(bam_to_h5_cmd_template)
  print(cmd)
  if (FALSE) # TODO uncomment
  {
  exit_code <- system(cmd)
  print(glue("bam_to_h5.R exit code: {exit_code}"))
  expect_equal(exit_code, 0, info = "Unexpected exit code from bam_to_h5.R")
  }

  ##### EXTRACT GFF (generic) #####

  ## bam_to_h5.R-style

  print("========== GFF - bam_to_h5.R-style ==========")
  gff <- rtracklayer::readGFFAsGRanges(gff_file)
  print(gff)
  print(length(gff)) # 204
  print(class(gff)) # GRanges attr(,"package") GenomicRanges
  print(typeof(gff)) # S4
  gff_names <- mcols(gff)
  print(gff_names)
  print(length(gff_names)) # 5 
  print(class(gff_names)) # DFrame attr(,"package") S4Vectors
  print(typeof(gff_names)) # S4
  print(dim(gff_names)) # 204 5
  gff_names <- gff_names["Name"] # --primary_id  
  print(gff_names)
  gff_names <- gff_names[,1] 
  print(gff_names)
  gff_names <- unique(gff_names)
  print(gff_names)
  print(length(gff_names)) # 68
  print(class(gff_names)) # character
  print(typeof(gff_names)) # character
  gff_pid <- mcols(gff)["Name"][,1]
  print(gff_pid)
  print(length(gff_pid)) # 204
  print(class(gff_pid)) # character
  print(typeof(gff_pid)) # character

  # gene <- "YAL062W"
  # gene_location <- gff[gff_pid == gene]
  # print(gene_location)
  # gene_location <- gff["Name" == gene]
  # print(gene_location)

  ## read_count_functions.R-style

  print("========== GFF ==========")
  gff_df <- readGFFAsDf(gff_file)
  print(gff_df)
  gff_names <- unique(gff_df$seqnames) # Tibble data frame
  print(gff_names)
  # [1] YAL068C YAL067W-A ...
  # [64] ... YAL002W YAL001C  
  # 68 Levels: YAL001C YAL002W YAL003W YAL005C YAL007C YAL008W YAL009W ... YAL068C
  print(length(gff_names)) # 68
  print(class(gff_names)) # factor
  print(typeof(gff_names)) # integer
  print(levels(gff_names))
  # [1] "YAL001C"   "YAL002W"  ...
  # [67] "YAL067W-A" "YAL068C"  
 
  ##### EXTRACT BAM (generic) #####

  print("========== BAM ==========")

  # https://kasperdanielhansen.github.io/genbioconductor/html/Rsamtools.html
  bam_file_f <- BamFile(bam_file)
  print(bam_file_f)
  bam_hdr_seq_info <- seqinfo(bam_file_f)
  print(bam_hdr_seq_info)
  #Seqinfo object with 68 sequences from an unspecified genome:
  #seqnames  seqlengths isCircular genome
  #YAL068C          863       <NA>   <NA>
  # ...
  #YAL001C         3983       <NA>   <NA>
  print(class(bam_hdr_seq_info)) # SeqInfo. attr(,"package"), GenomeInfoDb
  print(typeof(bam_hdr_seq_info)) # S4
  print(length(bam_hdr_seq_info)) # 68
  print(countBam(bam_file_f))
  #  space start end width       file records nucleotides
  # 1    NA    NA  NA    NA WTnone.bam   14516      399027

  bam_hdr_seq_names <- bam_hdr_seq_info@seqnames
  print(bam_hdr_seq_names)
  print(length(bam_hdr_seq_names)) # 68
  print(class(bam_hdr_seq_names)) # character
  print(typeof(bam_hdr_seq_names)) # character
  # [1] "YAL068C"   "YAL067W-A"  ...
  # [67] "YAL002W" "YAL001C"  

  # Alternative approach
  # https://www.rdocumentation.org/packages/GenomicAlignments/versions/1.8.4/topics/readGAlignments
  # Uses bam <- readGAlignments(bam_file)
  # By default readGAlignments extracts:
  # seqnames, strand, cigar, qwidth, start, end, width
  # We also want "flag" so need to specify it explicitly.
  bam_what <- c("flag")
  bam_params <- ScanBamParam(what=bam_what)
  bam <- readGAlignments(bam_file, param=bam_params, use.names=T)
  print(class(bam)) # GAlignments attr(,"package") GenomicAlignments
  print(typeof(bam)) # S4
  print(length(bam)) # 14516
  print(bam)
  # GAlignments object with 14516 alignments and 1 metadata column:
  #          seqnames strand       cigar    qwidth     start       end     width
  #             <Rle>  <Rle> <character> <integer> <integer> <integer> <integer>
  #      [1]  YAL062W      +       26M2S        28       752       777        26
  #      [2]  YAL062W      +       25M2S        27       753       777        25
  #      ...       ... .       ...
  #  [14512]  YAL001C      +         26M        26       539       564        26
  #  [14513]  YAL001C      +         15M        15      1753      1767        15
  #  [14514]  YAL001C      +       1S27M        28      3007      3033        27
  #  [14515]  YAL001C      -       28M1S        29      3559      3586        28
  #  [14516]  YAL001C      -       1S12M        13      3562      3573        12
  #              njunc |      flag
  #          <integer> | <integer>
  #      [1]         0 |         0
  #      [2]         0 |         0
  #      ...       ... .       ...
  #  [14512]         0 |         0
  #  [14513]         0 |         0
  #  [14514]         0 |         0
  #  [14515]         0 |       272
  #  [14516]         0 |        16
  print(seqnames(bam))
  print(length(seqnames(bam))) # 14516
  print(length(sort(seqnames(bam)))) # 14516
  bam_seq_names <- unique(sort(seqnames(bam)))
  print(class(bam_seq_names)) #  factor
  print(typeof(bam_seq_names)) # integer
  print(length(bam_seq_names)) # 57

  ##### EXTRACT AND VALIDATE H5 (generic) #####

  print("========== H5 ==========")

  h5_data <- rhdf5::h5ls(h5_file, recursive = 1)
  print(h5_data)
  #  group      name     otype dclass dim
  # 0      /   YAL001C H5I_GROUP           
  # 58     /   YAL062W H5I_GROUP
  print(length(h5_data)) # 5
  print(class(h5_data)) # data.frame
  print(typeof(h5_data)) # list
  print(dim(h5_data)) # 68 5
  h5_names <- h5_data$name
  print(h5_names)
  # [1] "YAL001C"   "YAL002W"  ...
  # [67] "YAL067W-A" "YAL068C"  
  print(length(h5_names)) # 68
  print(class(h5_names)) # character
  print(typeof(h5_names)) # character

  # Validate against GFF
  expect_equal(length(h5_names), length(gff_names),
    info = "Unexpected number of sequence names when compared to GFF")
  expect_equal(as.factor(sort(h5_names)), sort(gff_names),
    info = "Unexpected sequence names when compared to GFF")

  # Validate against BAM header
  expect_equal(length(h5_names), length(bam_hdr_seq_names),
    info = "Unexpected number of sequence names when compared to BAM header")
  expect_equal(sort(h5_names), sort(bam_hdr_seq_names),
    info = "Unexpected sequence names when compared to BAM header")

  # Validate againt BAM content
  expect_true(length(bam_hdr_seq_names) <= length(h5_names),
    info = "Expected number of sequence names to be greater or equal to those in BAM")
  expect_true(all(sort(bam_seq_names) %in% as.factor(sort(h5_names))),
    info = "Expected sequence names superset of those in BAM")

  ##### EXTRACT GFF (gene-specific) #####

  print("Genes")
  for (gene in h5_names)
  {
    print(gene)
  }
  gene <- "YAL062W"
  # gene <- "YAL001C"
  print(gene)
  print(gff_df)
  print(class(gff_df)) # tbl_df tbl data.frame
  print(typeof(gff_df)) # list
  print(dim(gff_df)) # 204 10
  gff_utr5_start <- GetCDS5start(gene, gff_df, ftype="UTR5")
  gff_utr5_end <- GetCDS3end(gene, gff_df, ftype="UTR5")
  gff_utr5_length <- gff_utr5_end - gff_utr5_start + 1
  gff_cds_start <- GetCDS5start(gene, gff_df, ftype="CDS")
  gff_cds_end <- GetCDS3end(gene, gff_df, ftype="CDS")
  gff_cds_length <- gff_cds_end - gff_cds_start + 1
  gff_utr3_start <- GetCDS5start(gene, gff_df, ftype="UTR3")
  gff_utr3_end <- GetCDS3end(gene, gff_df, ftype="UTR3")
  gff_utr3_length <- gff_utr3_end - gff_utr3_start + 1
  print("UTR5 start/end:")
  print(gff_utr5_start) # 1
  print(gff_utr5_end) # 250
  print(gff_utr5_length) # 250
  print("CDS start/end:")
  print(gff_cds_start) # 251
  print(gff_cds_end) # 1624 
  print(gff_cds_length) # 1374
  print("UTR3 start/end:")
  print(gff_utr3_start) # 1625
  print(gff_utr3_end) # 1874
  print(gff_utr3_length) # 250

  ##### EXTRACT BAM (gene-specific) #####

  bam_hdr_gene <- bam_hdr_seq_info[gene]
  print(bam_hdr_gene)
  # seqnames seqlengths isCircular genome
  # YAL062W        1874         NA   <NA>
  print(class(bam_hdr_gene)) # SeqInfo. attr(,"package"), GenomeInfoDb
  print(typeof(bam_hdr_gene)) # S4
  bam_hdr_gene_seq_length <- bam_hdr_gene@seqlengths
  print(bam_hdr_gene_seq_length) # 1874

  # https://www.rdocumentation.org/packages/GenomicAlignments/versions/1.8.4/topics/readGAlignments
  bam_gene = bam[(seqnames(bam) == gene)]
  bam_gene = bam[(seqnames(bam) == gene)]
  print(bam_gene)
  # GAlignments object with 2 alignments and 1 metadata column:
  #       seqnames strand       cigar    qwidth     start       end     width
  #          <Rle>  <Rle> <character> <integer> <integer> <integer> <integer>
  #   [1]  YAL062W      +       26M2S        28       752       777        26
  #   [2]  YAL062W      +       25M2S        27       753       777        25
  #           njunc |      flag
  #       <integer> | <integer>
  #   [1]         0 |         0
  #   [2]         0 |         0
  #   seqinfo: 68 sequences from an unspecified genome
  # For YAL001C:
  #   SRR1042855.35576357  YAL001C      +       2S28M        30       162       189
  #   SRR1042855.43554901  YAL001C      +         26M        26       539       564
  #   SRR1042855.18823368  YAL001C      +         15M        15      1753      1767
  #   SRR1042855.38021801  YAL001C      +       1S27M        28      3007      3033
  #   SRR1042855.35349348  YAL001C      -       28M1S        29      3559      3586
  #   SRR1042855.43963789  YAL001C      -       1S12M        13      3562      3573
  #   SRR1042855.35576357        28         0 |         0
  #   SRR1042855.43554901        26         0 |         0
  #   SRR1042855.18823368        15         0 |         0
  #   SRR1042855.38021801        27         0 |         0
  #   SRR1042855.35349348        28         0 |       272
  #   SRR1042855.43963789        12         0 |        16
  print(names(bam_gene)) # "SRR1042855.5473767" "SRR1042855.1850623"
  print(class(bam_gene)) # GAlignments attr(,"package") GenomicAlignments
  print(typeof(bam_gene)) # S4
  print(length(bam_gene)) # 2. For YAL001C 6
  print(mcols(bam_gene))
  #  DataFrame with 2 rows and 1 column
  #                        flag
  #                   <integer>
  # SRR1042855.5473767         0
  # SRR1042855.1850623         0
  print(mcols(bam_gene)$flag) # 0 0. For YAL001C 0 0 0 0 272 16
  print("Sequences with flag == 0")
  bam_gene_flag_zero = bam[(seqnames(bam) == gene) & (mcols(bam)$flag == 0)]
  print(bam_gene_flag_zero)
  # GAlignments object with 2 alignments and 1 metadata column:
  #       seqnames strand       cigar    qwidth     start       end     width
  #          <Rle>  <Rle> <character> <integer> <integer> <integer> <integer>
  #   [1]  YAL062W      +       26M2S        28       752       777        26
  #   [2]  YAL062W      +       25M2S        27       753       777        25
  #           njunc |      flag
  #       <integer> | <integer>
  #   [1]         0 |         0
  #   [2]         0 |         0
  #   seqinfo: 68 sequences from an unspecified genome
  # For YAL001C:
  #   SRR1042855.35576357  YAL001C      +       2S28M        30       162       189
  #   SRR1042855.43554901  YAL001C      +         26M        26       539       564
  #   SRR1042855.18823368  YAL001C      +         15M        15      1753      1767
  #   SRR1042855.38021801  YAL001C      +       1S27M        28      3007      3033
  #   SRR1042855.35576357        28         0 |         0
  #   SRR1042855.43554901        26         0 |         0
  #   SRR1042855.18823368        15         0 |         0
  #   SRR1042855.38021801        27         0 |         0
  print(length(bam_gene_flag_zero)) # 2. For YAL001C 4
  print("Sequences with flag != 0")
  bam_gene_flag_non_zero = bam[(seqnames(bam) == gene) & (mcols(bam)$flag != 0)]
  print(bam_gene_flag_non_zero)
  # GAlignments object with 0 alignments and 1 metadata column:
  #    seqnames strand       cigar    qwidth     start       end     width
  #       <Rle>  <Rle> <character> <integer> <integer> <integer> <integer>
  #        njunc |      flag
  #    <integer> | <integer>
  # For YAL001C:
  #   SRR1042855.35349348  YAL001C      -       28M1S        29      3559      3586
  #   SRR1042855.43963789  YAL001C      -       1S12M        13      3562      3573
  #   SRR1042855.35349348        28         0 |       272
  #   SRR1042855.43963789        12         0 |        16
  print(length(bam_gene_flag_non_zero)) # 0. For YAL001C 2

  ##### EXTRACT AND VALIDATE H5 (gene-specific) #####

  # 'buffer_left': number of nucleotides upstream of the start codon (ATG) (UTR5 length) (from bam_to_h5.R command-line)
  print("buffer_left:")
  h5_buffer_left <- GetGeneBufferLeft(gene, dataset, h5_file) # double
  print(h5_buffer_left) # 250
  expect_equal(h5_buffer_left, buffer,
    info = "Unexpected buffer_left when compared to bam_to_h5.R parameter")
  expect_equal(h5_buffer_left, gff_utr5_length,
    info = "Unexpected buffer_left when compared to GFF UTR5 length")

  # 'buffer_right': number of nucleotides downstream of the stop codon (TAA/TAG/TGA) (UTR3 length) (from bam_to_h5.R command-line)
  print("buffer_right:")
  h5_buffer_right <- GetGeneBufferRight(gene, dataset, h5_file) # integer
  print(h5_buffer_right) # 250
  expect_equal(h5_buffer_right, buffer,
    info = "Unexpected buffer_right when compared to bam_to_h5.R parameter")
  expect_equal(h5_buffer_left, gff_utr3_length,
    info = "Unexpected buffer_left when compared to GFF UTR3 length")

  # 'start_codon_pos': Positions corresponding to start codon of CDS in organism sequence (from GFF)
  expected_start_codons <- as.array(seq(gff_cds_start, gff_cds_start + 2))
  print("expected_start_codons:")
  print(expected_start_codons) # 251, 252, 253
  print("start_codon_pos:")
  h5_start_codon_pos <- GetGeneStartCodonPos(gene, dataset, h5_file) # 1D array of 3 integer
  print(h5_start_codon_pos) # 251 252 253
  expect_equal(length(h5_start_codon_pos), 3,
    info = "Unexpected number of start_codon_pos, expected 3")
  expect_equal(h5_start_codon_pos, expected_start_codons,
    info = "Unexpected start_codon_pos when compared to GFF CDS start codon positions")

  # 'stop_codon_pos': Positions corresponding to stop codon of CDS in organism sequence (from GFF)
  expected_stop_codons <- as.array(seq(gff_cds_end - 2, gff_cds_end))
  print("expected_stop_codons:")
  print(expected_stop_codons) # 1622 1623 1624
  print("stop_codon_pos:")
  h5_stop_codon_pos <- GetGeneStopCodonPos(gene, dataset, h5_file) # 1D array of 3 integer
  print(h5_stop_codon_pos) # 1622 1623 1624
  expect_equal(length(h5_stop_codon_pos), 3,
    info = "Unexpected number of stop_codon_pos, expected 3")
  expect_equal(h5_stop_codon_pos, expected_stop_codons,
    info = "Unexpected stop_codon_pos when compared to GFF CDS stop codon positions")

  # 'lengths' : Lengths of mapped reads.
  print("lengths:")
  expected_lengths <- as.array(seq(min_read_length, max_read_length))
  h5_lengths <- GetGeneMappedReadLengths(gene, dataset, h5_file) # 1D array of <max_read_length - min_read_length + 1> integer
  print(h5_lengths)
  # [1] 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
  # [26] 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
  expect_equal(length(h5_lengths), num_read_counts,
    info = "Number of lengths does not equal max_read_length - min_read_length + 1")
  expect_equal(h5_lengths, expected_lengths,
      info = "Unexpected lengths")

  # 'reads_by_len': Counts of number of ribosome sequences of each length (from BAM).
  print("reads_by_len:")
  h5_reads_by_len <- GetGeneReadLength(gene, dataset, h5_file) # 1D array of <max_read_length - min_read_length + 1> double
  print(h5_reads_by_len)
  # [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  # [39] 0 0 0
  expect_equal(length(h5_reads_by_len), num_read_counts,
    info = "Number of reads_by_len does not equal max_read_length - min_read_length + 1")
  # Calculate expected reads_by_len based on information from BAM
  expected_reads_by_len <- as.array(replicate(num_read_counts, 0))
  print(expected_reads_by_len)
  # [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  # [39] 0 0 0
  print(length(expected_reads_by_len)) # 41
  for (width in sort(qwidth(bam_gene_flag_zero)))
  {
      width_index <- width - min_read_length + 1
      expected_reads_by_len[width_index] <- expected_reads_by_len[width_index] + 1
  }
  print(expected_reads_by_len)
  # [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  # [39] 0 0 0
  expect_equal(h5_reads_by_len, expected_reads_by_len,
    info = "reads_by_len does not correspond with widths in BAM")

  # 'reads_total': Total number of ribosome sequences (from BAM, equal to number of non-zero reads in 'reads_by_len').
  h5_reads_total <- GetGeneReadsTotal(gene, dataset, h5_file) # 1D array of 1 double
  print("reads_total:")
  print(h5_reads_total) # 2
  h5_reads_len_total <- Reduce("+", h5_reads_by_len)
  print("SUM(reads_by_len):")
  print(h5_reads_len_total)
  expect_equal(h5_reads_total[1], h5_reads_len_total,
    info = "reads_total does not equal sum of totals in reads_by_len")
  expected_reads_total = length(bam_gene_flag_zero)
  print("Number of Flag = 0 BAM records:")
  print(expected_reads_total) # 2
  expect_equal(length(h5_reads_total), 1,
      info = "Unexpected number of reads_total")
  expect_equal(h5_reads_total[1], expected_reads_total,
      info = "reads_total does not equal number of BAM records with Flag = 0")

  # 'data': Positions and lengths of ribosome sequences within the organism data (from BAM).
  print("data:")
  h5_data <- GetGeneDatamatrix(gene, dataset, h5_file)
  # print(h5_data) # Verbose
  print(length(h5_data)) # 76384
  print(class(h5_data)) # matrix
  print(typeof(h5_data)) # integer
  print(dim(h5_data)) # 41 1874
  print(nrow(h5_data)) # 41
  print(ncol(h5_data)) # 1874
  expect_equal(nrow(h5_data), num_read_counts,
    info = "Number of data rows does not equal max_read_length - min_read_length + 1")
  expected_num_data_cols <- h5_stop_codon_pos[3] + buffer
  print(expected_num_data_cols) # 1874
  expect_equal(ncol(h5_data), expected_num_data_cols,
    info = "Number of data columns does not equal stop_codon_pos[3] + buffer")
  expect_equal(ncol(h5_data), gff_utr3_end,
    info = "Number of data columns does not equal GFF UTR3 final nt position")
  expect_equal(ncol(h5_data), bam_hdr_gene_seq_length,
    info = "Number of data columns does not equal sequence lengths from BAM header")

  # TODO Generalise - see below.
  expected_row <- integer(num_read_counts)
  expected_row[19] <- 1 #  1-indexed R vs 0-indexed H5 (18)
  row <- h5_data[,752] # 1-indexed R vs 0-indexed H5
  print(row)
  expect_equal(row, expected_row, info = "Unexpected data[751]")

  expected_row <- integer(num_read_counts)
  expected_row[18] <- 1 #  1-indexed R vs 0-indexed H5 (17)
  row <- h5_data[,753] # 1-indexed R vs 0-indexed H5
  print(row)
  expect_equal(row, expected_row, info = "Unexpected data[752]")

  ##### EXTRACT BAM WIP #####
  ##### EXTRACT BAM WIP #####
  ##### EXTRACT BAM WIP #####

  # TODO
  # TODO
  # TODO
  
  print(gene)
  print(bam_hdr_gene)
  print(bam_gene_flag_zero)
  #                     seqnames strand       cigar    qwidth     start       end
  #                        <Rle>  <Rle> <character> <integer> <integer> <integer>
  #  SRR1042855.5473767  YAL062W      +       26M2S        28       752       777
  #  SRR1042855.1850623  YAL062W      +       25M2S        27       753       777
  #                         width     njunc |      flag
  #                     <integer> <integer> | <integer>
  #  SRR1042855.5473767        26         0 |         0
  #  SRR1042855.1850623        25         0 |         0

  # TODO Use bam_gene_flag_zero, start, qwidth

  # SRR1042855.5473767	0	YAL062W	752	60	26M2S	*	0	0	TCATACAAGAACTCCTGGGAAGGTGTCT	CCCFFFFFHHHHHJJJJJJJJJJEGHIJ	AS:i:-2XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:26	YT:Z:UU	XS:A:+	NH:i:1
  # Position 752 Length 28
  # SRR1042855.1850623	0	YAL062W	753	60	25M2S	*	0	0CATACAAGAACTCCTGGGAAGGTGTCT	@@@DDDDDHHH?FGFIIDHGII+ACDH	AS:i:-2	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:25	YT:Z:UU	XS:A:+	NH:i:1
 # Position 753 Length 27

  # 'data': Positions and lengths of ribosome sequences within the organism data (from BAM).
  # TODO create expected "data"
  # TODO reads_by_len[i] sum of DATA[*, i] sum over all positions for a specific length.
  # TODO Check DATA[p, i] = 1 if there is a sequence from BAM at position p+1 which has length equal to lengths[i], else 0.
  # TODO check sequence with "non-zeros" is in BAM.
  # TODO check sequence with "zeros" only is not in BAM.
  # TODO Cross-check reads_by_len[i] = sum of DATA[*, i] i.e. sum across all positions for a specific length.

  # TODO Try YAL001C
  # TODO Extend to iterate through all sequences in H5.

  expect_equal(0, 0, info = "Example assertion") # TODO remove
})
