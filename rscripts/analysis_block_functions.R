# analysis_block_functions.R
# R functions riboviz analysis blocks
# https://github.com/riboviz/riboviz/


#####

#
#
# START 3NT PERIODICITY
#
#

# MEDIUM FUNCTIONS:

CalculateThreeNucleotidePeriodicity <- function(gene_names, dataset, hd_file, gff_df){
  
  # get gene and position specific total counts for all read lengths
  gene_poslen_counts_5start_df <- AllGenes5StartPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hd_file = hd_file, gff_df = gff_df)
  
  gene_poslen_counts_3end_df <- AllGenes3EndPositionLengthCountsTibble(gene_names = gene_names, dataset= dataset, hd_file = hd_file, gff_df = gff_df)
  
  # summarize by adding different read lengths
  gene_pos_counts_5start <- gene_poslen_counts_5start_df %>%
    group_by(Pos) %>%
    summarize(Counts = sum(Counts))
  # gives:
  # > str(gene_pos_counts_5start)
  # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	75 obs. of  2 variables:
  #   $ Pos   : int  -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 ...
  #   $ Counts: int  285 318 307 386 291 347 840 330 475 355 ...
  
  gene_pos_counts_3end <- gene_poslen_counts_3end_df  %>%
    group_by(Pos) %>%
    summarize(Counts = sum(Counts))
  # gives:
  # > str(gene_pos_counts_3end)
  # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	75 obs. of  2 variables:
  #   $ Pos   : int  -49 -48 -47 -46 -45 -44 -43 -42 -41 -40 ...
  #   $ Counts: int  19030 13023 50280 19458 12573 46012 19043 13282 36968 20053 ...
  
  three_nucleotide_periodicity_data <- bind_rows(
    gene_pos_counts_5start %>% mutate(End = "5'"),
    gene_pos_counts_3end %>% mutate(End = "3'")
  ) %>%
    mutate(End = factor(End, levels = c("5'", "3'")))
  # gives:
  # > str(three_nucleotide_periodicity_data)
  # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	150 obs. of  3 variables:
  #   $ Pos   : int  -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 ...
  #   $ Counts: int  285 318 307 386 291 347 840 330 475 355 ...
  #   $ End   : Factor w/ 2 levels "5'","3'": 1 1 1 1 1 1 1 1 1 1 ...
  
  return(three_nucleotide_periodicity_data)
  
} # end CalculateThreeNucleotidePeriodicity() definition
# gives:
#   CalculateThreeNucleotidePeriodicity(gene_names = gene_names, dataset = dataset, hd_file = hd_file, gff_df = gff_df)
#   # A tibble: 150 x 3
#   Pos Counts End
#   <int>  <int> <fct>
#     1   -24    285 5'
#     2   -23    318 5'
#     3   -22    307 5'
#     4   -21    386 5'
#     5   -20    291 5'
#     6   -19    347 5'
#     7   -18    840 5'
#     8   -17    330 5'
#     9   -16    475 5'
#    10   -15    355 5'
#   # … with 140 more rows

# define PlotThreeNucleotidePeriodicity() function with reasonable arguments
PlotThreeNucleotidePeriodicity <- function(three_nucleotide_periodicity_data){
  
  # Plot
  three_nucleotide_periodicity_plot <- ggplot(
    three_nucleotide_periodicity_data,
    aes(x = Pos, y = Counts)) +
    geom_line() +
    facet_wrap(~End, scales = "free") +
    labs(x = "Nucleotide Position", y = "Read counts")
  
  return(three_nucleotide_periodicity_plot)
  
} # end PlotThreeNucleotidePeriodicity() definition

# potentially replace/tweak plot_ribogrid() to follow StyleGuide
PlotStartCodonRiboGrid <- function(gene_poslen_counts_5start_df){
  # function to do the ribogrid & ribogridbar plots?
  # ribogrid_5start
  start_codon_ribogrid_plot <- plot_ribogrid(gene_poslen_counts_5start_df)
  return(start_codon_ribogrid_plot)
} # end PlotStartCodonRiboGrid() definition

SaveStartCodonRiboGrid <- function(start_codon_ribogrid_plot){
  # function to do the ribogrid & ribogridbar plots?
  # ribogrid_5start
  start_codon_ribogrid_plot %>%
    ggsave(
      filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogrid.pdf")),
      width = 6, height = 3
    )
  #return() # no return as writing-out
} # end SaveStartCodonRiboGrid() definition

PlotStartCodonRiboGridBar <- function(gene_poslen_counts_5start_df){
  start_codon_ribogrid_bar_plot <- barplot_ribogrid(gene_poslen_counts_5start_df)
  return(start_codon_ribogrid_bar_plot)
} # end PlotStartCodonRiboGridBar() definition

SaveStartCodonRiboGridBar <- function(start_codon_ribogrid_bar_plot){
  start_codon_ribogrid_bar_plot %>%
    ggsave(
      filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogridbar.pdf")),
      width = 6, height = 5
    )
  #return() # no return as writing-out
} # end SaveStartCodonRiboGridBar() definition

SavePlotThreeNucleotidePeriodicity <- function(three_nucleotide_periodicity_plot) {
  # Save plot and file
  ggsave(
    three_nucleotide_periodicity_plot,
    filename = file.path(output_dir, paste0(output_prefix, "3nt_periodicity.pdf"))
  )
  # return() # NO RETURN as writing out
} # end of function definition SavePlotThreeNucleotidePeriodicity()

WriteThreeNucleotidePeriodicity <- function(three_nucleotide_periodicity_data) {
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "3nt_periodicity.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    three_nucleotide_periodicity_data,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F)
  # return()? NO RETURN
}  # end of function definition WriteThreeNucleotidePeriodicity()


#
#
# END 3NT PERIODICITY
#
# START ALL MAPPED READS
#
#

# MEDIUM FUNCTIONS:

# calculate function
CalculateReadLengths <- function(gene_names, dataset, hd_file){
  
  ## distribution of lengths of all mapped reads
  print("Starting: Distribution of lengths of all mapped reads")
  
  # read length-specific read counts stored as attributes of 'reads' in H5 file
  gene_sp_read_length <- lapply(gene_names, function(gene) {
    GetGeneReadLength(gene, hd_file)
  })
  
  # sum reads of each length across all genes
  read_length_data <- data.frame(
    Length = read_range,
    Counts = gene_sp_read_length %>%
      Reduce("+", .)
  )
  
  # return read length data
  return(read_length_data)
  
} # end definition of function CalculateReadLengths()
# > str(read_length_data)
# 'data.frame':	41 obs. of  2 variables:
#   $ Length: int  10 11 12 13 14 15 16 17 18 19 ...
#   $ Counts: num  0 0 0 0 0 ...

# plot function
PlotReadLengths <- function(read_length_data){
  
  # plot read lengths with counts
  read_len_plot <- ggplot(read_length_data, aes(x = Length, y = Counts)) +
    geom_bar(stat = "identity")
  
  return(read_len_plot)
  
} # end definition of function PlotReadLengths()

# save pdf
SavePlotReadLengths <- function(read_len_plot) {
  
  ggsave(read_len_plot, filename = file.path(output_dir, paste0(output_prefix, "read_lengths.pdf")))
  
  # return() # NO RETURN as writing out
  
} # end definition of function SavePlotReadLenths()

# save read lengths plot and file
WriteReadLengths <- function(read_length_data){
  
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "read_lengths.tsv"))
  
  write_provenance_header(path_to_this_script, tsv_file_path)
  
  write.table(
    read_length_data,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
  # return() NO RETURN as writing out
} # end definition of function WriteReadLengths()


#
#
# END ALL MAPPED READS
#
# START BIASES IN NUCLEOTIDE COMPOSITION
#
#

# MEDIUM FUNCTIONS:


# CALCULATE Biases In Nucleotide Composition Along Mapped Read Lengths
CalculateBiasesInNucleotideComposition <- function(gene_names, dataset, hd_file, read_range, min_read_length){
  
  # This is in a conditional loop because it fails for some inputs
  # and has not been debugged. Needs to be rewritten in tidyverse
  all_out <- c() # creates output object set to null
  for (lid in seq_len(length(read_range))) { # TODO: WHAT IS LID
    out <- lapply(gene_names, function(x) { # TODO: fix variablename 'out' (also check below for uses)
      # For each read length convert reads to IRanges
      GetNTReadPosition(gene = as.character(x),
                        dataset = dataset,
                        hd_file = hd_file,
                        lid = lid, min_read_length = min_read_length)
    })
    names(out) <- gene_names
    
    # read_range is sequence 10:50 (inclusive)
    # length(read_range) is 41
    # seq_len(41) generates sequence 1:41
    
    # GetNTReadPosition:
    # GetNTReadPosition <- function(gene, dataset, hd_file, lid, min_read_length) {
    #   reads_pos_len <- GetGeneDatamatrix(gene, dataset, hd_file)[lid, ] # Get reads of a particular length
    #   reads_pos_len <- reads_pos_len[1:(length(reads_pos_len) - (lid + min_read_length - 1))] # Ignore reads whose 5' ends map close to the end of the 3' buffer
    #   pos <- rep(1:length(reads_pos_len), reads_pos_len) # nt positions weighted by number of reads mapping to it
    #   pos_IR <- IRanges::IRanges(start = pos, width = (lid + min_read_length - 1)) # Create an IRanges object for position-specific reads of a particular length
    #   return(pos_IR)
    # }
    
    # TODO FLIC testcode
    # for (lid in 1:3) { # TODO: WHAT IS LID
    #   out <- lapply(gene_names, function(x) { # TODO: fix variablename 'out' (also check below for uses)
    #     # For each read length convert reads to IRanges
    #     GetNTReadPosition(gene = as.character(x),
    #                       dataset = dataset,
    #                       hd_file = hd_file,
    #                       lid = lid, min_read_length = min_read_length)
    #   })
    # names(out) <- gene_names[1]
    
    # Get position-specific nucleotide counts for reads in each frame
    
    # # TODO: test this function
    # CalculatePositionSpecificNucleotideCountsByFrame <- function(gene_names, cframe, lid, num_processes){
    #   mcapply(gene_names, function(gene){
    #     cons_mat(gene = gene, pos_IR = out[[gene]], cframe = cframe, lid = lid)
    #   }, mc.cores = num_processes)
    # } # end of function definition CalculatePositionSpecificNucleotideCountsByFrame()
    # test_frame0 <- CalculatePositionSpecificNucleotideCountsByFrame(gene_names, cframe=0, lid, num_processes)
    
    # frame 0
    fr0 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 0, lid = lid)
    }, mc.cores = num_processes)
    allfr0 <- do.call(rbind, fr0)
    
    # frame 1
    fr1 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 1, lid = lid)
    }, mc.cores = num_processes)
    allfr1 <- do.call(rbind, fr1)
    
    # frame 2
    fr2 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 2, lid = lid)
    }, mc.cores = num_processes)
    allfr2 <- do.call(rbind, fr2)
    
    # Get position-specific freq for all nucleotides
    cnt_fr0 <- signif(comb_freq(allfr0), 3)
    cnt_fr1 <- signif(comb_freq(allfr1), 3)
    cnt_fr2 <- signif(comb_freq(allfr2), 3)
    
    output <- data.frame(rbind(cnt_fr0, cnt_fr1, cnt_fr2)) # TODO: replace variablename 'output'
    all_out <- rbind(all_out, output)
  } # end of the for() loop
  
  # TODO:
  print("finished for loop")
  
  # Prepare variables for output file
  Length <- unlist(lapply(read_range, function(x) {
    rep(x, x * 3)
  }))
  Position <- unlist(lapply(read_range, function(x) {
    rep(1:x, 3)
  }))
  Frame <- unlist(lapply(read_range, function(x) {
    rep(0:2, each = x)
  }))
  
  # TODO:
  print("finished prepping Length, Position, Frame variables")
  
  all_out <- cbind(Length, Position, Frame, all_out)
  all_out[is.na(all_out)] <- 0
  
  # TODO:
  print("finished, returning all_out")
  
  return(all_out) # TODO: ensure this is correct; rename the returned output
} # end of function definition of CalculateBiasesInNucleotideComposition()
# run:
#CalculateBiasesInNucleotideComposition(gene_names, dataset, hd_file, read_range, min_read_length)

## WRITE DATA Biases In Nucleotide Composition Along Mapped Read Lengths
WriteBiasesInNucleotideComposition <- function(all_out){
  # save file
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "pos_sp_nt_freq.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(all_out, file = tsv_file_path, append = T, sep = "\t", row = F, col = T, quote = F)
} # end of function definition: WriteBiasesInNucleotideComposition()
# run:
#WriteBiasesInNucleotideComposition(all_out)

#
#
# END BIASES IN NUCLEOTIDE COMPOSITION
#
# START CALCULATE READ FRAME FOR EVERY ORF
#
#

## calculate read frame for every annotated ORF

# MEDIUM FUNCTIONS:

# MEDIUM FUNCTIONS:

# get asite displacement legnths
ReadAsiteDisplacementLengthFromFile <- function(asite_disp_length_file){
  asite_displacement_length <- readr::read_tsv(asite_disp_length_file,
                                               comment = "#"
  )
  return(asite_displacement_length)
}

CalculateGeneReadFrames <- function(dataset, hd_file, gff_df, min_read_length, asite_displacement_length_from_file) {
  # TODO: wrap in function
  gene_read_frames_data <- gff_df %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::select(gene = seqnames, left = start, right = end) %>%
    purrr::pmap_dfr(GetGeneReadFrame,
                    hd_file = hd_file,
                    dataset = dataset,
                    min_read_length = min_read_length,
                    asite_displacement_length = asite_displacement_length_from_file
    )
  return(gene_read_frames_data)
} # end CalculateGeneReadFrames() definition
# gives:
# TODO

PlotGeneReadFrames <- function(gene_read_frames_data){
  gene_read_frame_plot <- gene_read_frames_data %>%
    filter(Ct_fr0 + Ct_fr1 + Ct_fr2 > count_threshold) %>%
    BoxplotReadFrameProportion()
  
  return(gene_read_frame_plot)
} # end PlotGeneReadFrames() definition
# gives:
# TODO

SaveGeneReadFrames <- function(gene_read_frame_plot){
  # save read lengths plot and file
  ggsave(gene_read_frame_plot,
         filename = file.path(output_dir, paste0(output_prefix, "3ntframe_propbygene.pdf")),
         width = 3, height = 3
  )
  #return() # no return as writing-out
} # end SaveGeneReadFrames() definition
# gives:
# TODO

WriteGeneReadFrames <- function(gene_read_frames_data){
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "3ntframe_bygene.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    gene_read_frames_data,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
  #return() # no return as writing-out
} # end WriteGeneReadFrames() definition
# gives:
# TODO

#
#
# END CALCULATE READ FRAME FOR EVERY ORF
#
# START RPF POSITION SPECIFIC DISTRIBUTION OF READS
#
#

# codon-specific reads for RPF datasets
GetCodonPositionReads <- function(hd_file, gene, dataset, left, right, min_read_length) {
  # @ewallace: this needs documentation of inputs and outputs
  lid <- 28 - min_read_length + 1
  reads_pos <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
  reads_pos_subset <- reads_pos[, left:(dim(reads_pos)[2] - right)] # Subset positions such that only CDS codon-mapped reads are considered
  end_reads_pos_subset <- ncol(reads_pos_subset) # Number of columns of the subset
  
  l28 <- RcppRoll::roll_suml(reads_pos_subset[lid, 2:end_reads_pos_subset], n = 3, fill = NULL)[seq(1, length(reads_pos_subset[14, 2:end_reads_pos_subset]), 3)] # Map reads of length 28 to codons
  l29 <- RcppRoll::roll_suml(reads_pos_subset[(lid + 1), 2:end_reads_pos_subset], n = 3, fill = NULL)[seq(1, length(reads_pos_subset[15, 2:end_reads_pos_subset]), 3)] # Map reads of length 29 to codons
  l30 <- RcppRoll::roll_suml(reads_pos_subset[(lid + 2), 1:end_reads_pos_subset], n = 3, fill = NULL)[seq(1, length(reads_pos_subset[16, 1:end_reads_pos_subset]), 3)] # Map reads of length 30 to codons
  
  cod_sp_counts <- l28 + l29 + l30 # Sum of reads of lengths 28-30 at each codon
  cod_sp_counts <- cod_sp_counts[1:(length(cod_sp_counts) - 1)]
  return(cod_sp_counts)
}

# Nt-specific coverage for mRNA datasets
GetMRNACoverage <- function(hd_file, gene, dataset, left, right, read_range, min_read_length, buffer) {
  reads_pos <- GetGeneDatamatrix(gene, dataset, hd_file) # Get the matrix of read counts
  reads_pos_subset <- reads_pos[, left:(dim(reads_pos)[2] - right)] # Subset positions such that only CDS mapped reads are considered
  
  nt_IR_list <- lapply(read_range, function(w) {
    IRanges::IRanges(start = rep(1:ncol(reads_pos_subset), reads_pos_subset[(w - min_read_length + 1), ]), width = w)
  }) # Create list of IRanges for position-specific reads of all length
  nt_IR <- unlist(as(nt_IR_list, "IRangesList")) # Combine IRanges from different read lengths
  nt_cov <- IRanges::coverage(nt_IR) # Estimate nt-specific coverage of mRNA reads
  
  # Subset coverage to only CDS
  nt_counts <- rep.int(S4Vectors::runValue(nt_cov), S4Vectors::runLength(nt_cov))
  if (length(nt_counts) >= (buffer - left)) {
    nt_counts <- nt_counts[(buffer - left):length(nt_counts)]
  } else {
    nt_counts <- 0
  }
  
  cds_length <- ncol(reads_pos_subset) - (buffer - left - 1) # Length of CDS
  nt_sp_counts <- rep(0, cds_length)
  
  if (length(nt_counts) < cds_length) {
    if (length(nt_counts) > 0) {
      nt_sp_counts[1:length(nt_counts)] <- nt_counts
    }
  } else {
    nt_sp_counts <- nt_counts[1:cds_length]
  }
  return(nt_sp_counts)
}

#
#
# END RPF POSITION SPECIFIC DISTRIBUTION OF READS
#
# START MRNA POSITION SPECIFIC DISTRIBUTION OF READS
#
#


# # MEDIUM FUNCTIONS:



CalculateNucleotideBasedPositionSpecificReads <- function(gene, dataset, min_read_length, read_range, buffer){
  # create empty matrix to store position-specific read counts
  out5p <- matrix(NA, nrow = length(gene_names), ncol = 1500) # 5'
  out3p <- matrix(NA, nrow = length(gene_names), ncol = 1500) # 3'
  
  out <- lapply(gene_names, function(gene) {
    GetMRNACoverage(
      hd_file,
      gene,
      dataset,
      left = (buffer - 49),
      right = (buffer - 3),
      min_read_length = min_read_length,
      read_range,
      buffer
    )
  })
  names(out) <- gene_names
  
  cc <- 1
  for (gene in gene_names) {
    tmp <- out[[gene]]
    # only consider genes with at least 1 mapped read along its CDS
    if (sum(tmp) > 0) {
      tmp <- tmp / mean(tmp)
      if (length(tmp) > 1500) {
        out5p[cc, ] <- tmp[1:1500]
        out3p[cc, ] <- rev(tmp)[1:1500]
      } else {
        out5p[cc, 1:length(tmp)] <- tmp
        out3p[cc, 1:length(tmp)] <- rev(tmp)
      }
    }
    cc <- cc + 1
  }
  
  # estimate position-specific mean and std error of mapped read counts
  m5p <- signif(apply(out5p, 2, mean, na.rm = T), 4)
  m3p <- signif(apply(out3p, 2, mean, na.rm = T), 4)
  s5p <- signif(apply(out5p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)
  s3p <- signif(apply(out3p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)
  
  # normalize reads to last 150 nts of the 1500-nt window.
  # this allows easy comparison between datasets
  s5p <- s5p / mean(m5p[1350:1500])
  s3p <- s3p / mean(m3p[1350:1500])
  m5p <- m5p / mean(m5p[1350:1500])
  m3p <- m3p / mean(m3p[1350:1500])
  
  # create a dataframe to store the output for plots/analyses
  pos_sp_mrna_norm_coverage <- data.frame(
    Position = c(1:1500, 0:-1499),
    Mean = c(m5p, m3p),
    SD = c(s5p, s3p),
    End = factor(rep(c("5'", "3'"), each = 1500), levels = c("5'", "3'"))
  )
  
  return(pos_sp_mrna_norm_coverage)
  
} # end CalculateNucleotideBasedPositionSpecificReads() definition


PlotNucleotideBasedPositionSpecificReadsPerGene <- function(pos_sp_mrna_norm_coverage){
  # plot
  pos_sp_mrna_norm_coverage_plot <- ggplot(pos_sp_mrna_norm_coverage, aes(Position, Mean, col = End)) +
    geom_line() +
    facet_grid(~End, scales = "free") +
    guides(col = FALSE)
  
  return(pos_sp_mrna_norm_coverage_plot)
  
} # end PlotNucleotideBasedPositionSpecificReadsPerGene() definition

SaveNucleotideBasedPositionSpecificReadsPerGene <- function(pos_sp_mrna_norm_coverage_plot){
  # Save plot and file
  ggsave(pos_sp_mrna_norm_coverage_plot, filename = file.path(output_dir, paste0(output_prefix, "pos_sp_mrna_norm_coverage.pdf")))
  # return() # NO RETURN as writing out
  
}  # end SaveNucleotideBasedPositionSpecificReadsPerGene() definition

WriteNucleotideBasedPositionSpecificReadsPerGene <- function(pos_sp_mrna_norm_coverage){
  
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "pos_sp_mrns_norm_coverage.tsv"))
  write_provenance_header(this_script, tsv_file_path)
  write.table(
    pos_sp_mrna_norm_coverage,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
  # return() # NO RETURN as writing out
} # end WriteNucleotideBasedPositionSpecificReadsPerGene() definition

#
#
# END MRNA POSITION SPECIFIC DISTRIBUTION OF READS
#
# START TPMS OF GENES
#
#

# MEDIUM FUNCTIONS:

CalculateGeneTranscriptsPerMillion <- function(gene, dataset, hd_file){
  # calculate transcripts per million (TPM)
  # @FlicAnderson: what does reads_per_b stand for? Reads per Base? o.0  If so, maybe worth using full word?
  gene_sp_reads <- sapply(gene_names, GetGeneReadsTotal, dataset, hd_file)
  reads_per_b <- sapply(gene_names, GetGeneReadDensity, dataset, hd_file)
  
  tpms <- data.frame(
    ORF = gene_names,
    readcount = gene_sp_reads,
    rpb = reads_per_b,
    tpm = reads_per_b * 1e6 / sum(reads_per_b)
  )
  
  return(tpms)
} # end CalculateGeneTranscriptsPerMillion() definition


WriteGeneTranscriptsPerMillion <-  function(tpms){
  # write out to *_tpms.tsv
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "tpms.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    tpms,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
  # return() NO RETURN as writing out
} # end WriteGeneTranscriptsPerMillion() definition


#
#
# END TPMS OF GENES
#
# START TPMS CORRELATIONS WITH FEATURES
#
#

## Correlations between TPMs of genes with their sequence-based features

# MEDIUM FUNCTIONS:

# read features file
ReadSequenceBasedFeatures <- function(features_file){
  features <- read.table(features_file, h = T)
  return(features)
} # end ReadSequenceBasedFeatures() definition
# gives:
# TODO


CalculateSequenceBasedFeatures <- function(features, tpms){
  # Prepare data for plot
  # Consider only genes with at least count_threshold mapped reads
  features_plot_data <- merge(features, tpms, by = "ORF") %>%
    filter(readcount >= count_threshold, !is.na(ORF)) %>%
    select(-readcount, -rpb) %>%
    gather(Feature, Value, -ORF, -tpm)
  
  return(features_plot_data)
  
} # end CalculateSequenceBasedFeatures() definition
# gives:
# TODO


PlotSequenceBasedFeatures <- function(features_plot_data){
  
  features_plot <- ggplot(features_plot_data, aes(x = tpm, y = Value)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~Feature, scales = "free") +
    scale_x_log10() +
    geom_smooth(method = "lm") +
    xlab("TPM (transcripts per million)")
  
  return(features_plot)
} # end PlotSequenceBasedFeatures() definition


WriteSequenceBasedFeatures <- function(features_plot) {
  # Save plot and file
  ggsave(features_plot, filename = file.path(output_dir, paste0(output_prefix, "features.pdf")))
  
  # return() NO RETURN as writing out
  
} # end WriteSequenceBasedFeatures() definition


## Codon-specific ribosome densities for correlations with tRNAs


# # MEDIUM FUNCTIONS:
# 
# 
# CalculateCodonSpecificRibosomeDensity <- function(t_rna_file, codon_positions_file, gene, dataset, hd_file, buffer, count_threshold){
#   # dev version:
#   # yeast_tRNAs <- read.table(t_rna_file, h = T) # Read in yeast tRNA estimates
#   t_rna_df <- read.table(t_rna_file, h = T) # Read in yeast tRNA estimates
#   load(codon_positions_file) # Position of codons in each gene (numbering ignores first 200 codons)
#   # Reads in an object named "codon_pos"
#   out <- lapply(gene_names, function(gene) {
#     # From "Position specific distribution of reads" plot
#     GetCodonPositionReads(hd_file=hd_file, gene=gene, dataset=dataset,
#                           left = (buffer - 15), right = (buffer + 11),
#                           min_read_length = min_read_length)
#   }) # Get codon-based position-specific reads for each gene
#   names(out) <- gene_names
#   
#   gene_len <- sapply(out, length) # Calculate gene length in codons
#   out <- out[gene_len > 201] # Ignore genes with <=200 sense codons
#   
#   trim_out <- lapply(out, function(x) {
#     x[201:(length(x) - 1)]
#   }) # Trim first 200 codons and stop codon from each gene
#   read_counts_trim <- sapply(trim_out, sum) # Calculate read counts in trimmed genes
#   trim_out <- trim_out[read_counts_trim >= count_threshold] # Ignore genes with fewer than count_threshold mapped reads
#   
#   norm_out <- lapply(trim_out, function(x) {
#     x / mean(x)
#   }) # Normalize reads in each gene by their mean
#   
#   
#   # TODO: figure this out
#   # Calculate codon-specific mean ribosome-densities at A/P/E sites of the mapped reads
#   a_mn <- sapply(names(codon_pos), function(codon) {
#     mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
#       pos <- as.numeric(a[2])
#       norm_out[[a[1]]][pos]
#     })), na.rm = T)
#   })
#   p_mn <- sapply(names(codon_pos), function(codon) {
#     mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
#       pos <- as.numeric(a[2]) + 1
#       norm_out[[a[1]]][pos]
#     })), na.rm = T)
#   })
#   e_mn <- sapply(names(codon_pos), function(codon) {
#     mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
#       pos <- as.numeric(a[2]) + 2
#       norm_out[[a[1]]][pos]
#     })), na.rm = T)
#   })
#   
#   # Sort the values
#   A <- a_mn[order(names(codon_pos))]
#   P <- p_mn[order(names(codon_pos))]
#   E <- e_mn[order(names(codon_pos))]
#   
#   # TODO: this is a misnomer. Can calculate A/P/E-site norm density without t_rna_df.
#   # Should replace t_rna_df with an argument "codon_features_file"
#   # then plot against features in that, analogosly to (gene) to "features_file" above
#   cod_dens_tRNA_data <- cbind(t_rna_df, A, P, E)
#   
#   return(cod_dens_tRNA_data)
#   
# } # end of CalculateCodonSpecificRibosomeDensity() definition
# # gives:
# # TODO
# 

PlotCodonSpecificRibosomeDensityTRNACorrelation <- function(cod_dens_tRNA_data) {

  # Prepare data for plot
  cod_dens_tRNA_wide <- cod_dens_tRNA_data %>%
    gather(tRNA_type, tRNA_value, 3:6) %>%
    gather(Site, Ribodens, 3:5)

  # Plot
  cod_dens_tRNA_plot <- ggplot(cod_dens_tRNA_wide, aes(x = tRNA_value, y = Ribodens)) +
    geom_point(alpha = 0.3) +
    facet_grid(Site ~ tRNA_type, scales = "free_x") +
    geom_smooth(method = "lm") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(cod_dens_tRNA_plot)
} # end of PlotCodonSpecificRibosomeDensityTRNACorrelation() definition


SaveCodonSpecificRibosomeDensityTRNACorrelation <- function(cod_dens_tRNA_plot){

  # save plot
  ggsave(cod_dens_tRNA_plot, filename = file.path(output_dir, paste0(output_prefix, "codon_ribodens.pdf")))

  # return() # NO RETURN as writing out

} # end of SaveCodonSpecificRibosomeDensityTRNACorrelation() definition


WriteCodonSpecificRibosomeDensityTRNACorrelation <- function(cod_dens_tRNA_data){
  # Save file
  tsv_file_path <- file.path(output_dir, paste0(output_prefix, "codon_ribodens.tsv"))
  write_provenance_header(path_to_this_script, tsv_file_path)
  write.table(
    cod_dens_tRNA_data,
    file = tsv_file_path,
    append = T,
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
  # return() # NO RETURN as writing out
} # end of WriteCodonSpecificRibosomeDensityTRNACorrelation() definition

