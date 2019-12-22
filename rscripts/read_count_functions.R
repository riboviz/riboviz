# read_count_functions.R
# R functions for dealing with data in h5 "ribogrid" format
# https://github.com/riboviz/RiboViz/

# read in dependent packages
suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))
suppressMessages(library(rhdf5))
suppressMessages(library(parallel))
suppressMessages(library(optparse))
suppressMessages(library(RcppRoll))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(purrr))


#####
# Functions to read data from gff

readGFFAsDf <- purrr::compose(rtracklayer::readGFFAsGRanges, data.frame, as_tibble,
  .dir = "forward")

getCDS5start <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    filter(type==ftype, Name == name, strand == fstrand) %>% 
    pull(start) %>% 
    min 
}

getCDS3end <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    filter(type==ftype, Name == name, strand == fstrand) %>% 
    pull(end) %>% 
    max 
}

#####
# Functions to read data from h5 file

# function to get data matrix of read counts for gene and dataset from hdf5file
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

GetGeneReadLength <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Gopen(
      name = paste0("/", gene, "/", dataset, "/reads")
    ) %>%
      rhdf5::H5Aopen(name = "reads_by_len") %>%
    rhdf5::H5Aread() %>%
    return()
}

# function to get matrix of read counts from n_buffer before start codon to nnt_gene after
# for gene and dataset from hd5 file hdf5file, using UTR5 annotations in gff
GetGeneDatamatrix5start <- function(gene, dataset, hdf5file, 
                                    posn_5start,
                                    n_buffer, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  # @ewallace: replace this by gff_df?
  # n_utr5 <- BiocGenerics::width(gff[gff$type == "UTR5" & gff$Name == gene])
  # if n_buffer bigger than length n_utr5, pad with zeros:
  if (posn_5start > n_buffer) {
    # if posn_5start bigger than n_buffer
    n_left5 <- posn_5start - n_buffer # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  } else {
    # if length n_utr5 less than n_buffer
    n_left5 <- 1 # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = (n_buffer - posn_5start + 1 ))
  }
  n_right3 <- posn_5start + nnt_gene - 1 # column to end with (3'end)
  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  return(cbind(zeropad5_mat, data_mat_5start))
}

GetGeneDatamatrix3end <- function(gene, dataset, hdf5file, 
                                  posn_3end,
                                  n_buffer, nnt_gene) {
  # get data matrix of read counts from nnt_gene before stop codon to n_buffer after
  # for gene and dataset from hd5 file hdf5file, using UTR3 annotations in gff
  # if n_buffer bigger than length n_utr3, pad with zeros.
  # CHECK startpos/off-by-one
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_all <- ncol(data_mat_all)
  # @ewallace: replace this by gff_df?
  # n_utr3 <- BiocGenerics::width(gff[gff$type == "UTR3" & gff$Name == gene])
  n_left5 <- posn_3end - nnt_gene + 1 # column to start from (5'end)
  n_utr3  <- n_all - posn_3end
  if (n_utr3 >= n_buffer) {
    # length n_utr3 bigger than n_buffer
    n_right3 <- posn_3end + n_buffer # column to end with (3'end)
    zeropad3_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  } else {
    # length n_utr3 less than n_buffer
    n_right3 <- n_all # column to end with (3'end)
    zeropad3_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = n_buffer - n_utr3)
  }
  data_mat_3end <- data_mat_all[, n_left5:n_right3]
  return(cbind(data_mat_3end, zeropad3_mat))
}

TidyDatamatrix <- function(data_mat, startpos = 1, startlen = 1) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(data_mat) - 1)
  readlengths <- startlen:(startlen + nrow(data_mat) - 1)
  data_mat %>%
    set_colnames(positions) %>%
    as_tibble() %>%
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts))
}

plot_ribogrid <- function(tidymat) {
  ggplot(data = tidymat, aes(x = Pos, y = ReadLen, fill = Counts)) +
    geom_tile() +
    scale_fill_gradient("count", low = "white", high = "darkblue") +
    theme(panel.grid = element_blank()) +
    labs(x = "position of read 5' end", y = "read length")
}

barplot_ribogrid <- function(tidymat, small_read_range = 26:32) {
  ggplot(
    data = filter(tidymat, ReadLen %in% small_read_range),
    aes(x = Pos, y = Counts)
  ) +
    geom_col() +
    facet_grid(ReadLen ~ ., scales = "free_y") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "position of read 5' end", y = "count")
}

# GetGeneDatamatrix5start(gene="YAL003W",dataset="vignette",hdf5file,gff=gff) %>%
# TidyDatamatrix(startpos=-nnt_buffer,startlen=min_read_length) %>%
#   plot_ribogrid
# GetGeneDatamatrix3end(gene="YAL003W",dataset="vignette",hdf5file,gff=gff)

GetNTPeriod <- function(gene, dataset, hdf5file, left, right) {
  # previous version of script; not currently used
  data_mat <- GetGeneDatamatrix(gene, dataset, hdf5file)
  pos_sum <- colSums(data_mat) # Position-specific sum of reads of all lengths
  pos_sum <- pos_sum[left:(ncol(data_mat) - right)] # Ignore reads in parts of the buffer region defined by left/right
  return(pos_sum)
}

#####
## functions for position specific distribution of reads
# codon-specific reads for RPF datasets
GetCodonPositionReads <- function(gene, dataset, hdf5file, left, right, min_read_length) {
  # @ewallace: this needs documentation of inputs and outputs
  lid <- 28 - min_read_length + 1
  reads_pos <- GetGeneDatamatrix(gene, dataset, hdf5file) # Get the matrix of read counts
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
GetMRNACoverage <- function(gene, dataset, hdf5file, left, right, read_range, min_read_length, Buffer) {
  reads_pos <- GetGeneDatamatrix(gene, dataset, hdf5file) # Get the matrix of read counts
  reads_pos_subset <- reads_pos[, left:(dim(reads_pos)[2] - right)] # Subset positions such that only CDS mapped reads are considered

  nt_IR_list <- lapply(read_range, function(w) {
    IRanges::IRanges(start = rep(1:ncol(reads_pos_subset), reads_pos_subset[(w - min_read_length + 1), ]), width = w)
  }) # Create list of IRanges for position-specific reads of all length
  nt_IR <- unlist(as(nt_IR_list, "IRangesList")) # Combine IRanges from different read lengths
  nt_cov <- IRanges::coverage(nt_IR) # Estimate nt-specific coverage of mRNA reads

  # Subset coverage to only CDS
  nt_counts <- rep.int(S4Vectors::runValue(nt_cov), S4Vectors::runLength(nt_cov))
  if (length(nt_counts) >= (Buffer - left)) {
    nt_counts <- nt_counts[(Buffer - left):length(nt_counts)]
  } else {
    nt_counts <- 0
  }

  cds_length <- ncol(reads_pos_subset) - (Buffer - left - 1) # Length of CDS
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

#####
## functions to calculate read frame etc.

CalcAsiteFixedOneLength <- function(reads_pos_length, min_read_length,
                                    read_length, asite_disp) {
  # Calculate read A-site using a fixed displacement for a single read length
  length_row_choose <- read_length - min_read_length + 1
  reads_pos_length[length_row_choose, ] %>%
    dplyr::lag(n = asite_disp, default = 0)
}

CalcAsiteFixed <- function(reads_pos_length, min_read_length,
                           asite_disp_length = data.frame(
                             read_length = c(28, 29, 30),
                             asite_disp = c(15, 15, 15)
                           ),
                           colsum_out = TRUE) {
  # Calculate read A-site using a fixed displacement for fixed read lengths
  npos <- ncol(reads_pos_length)
  Asite_counts_bylength <-
    purrr::map2(
      asite_disp_length$read_length, asite_disp_length$asite_disp,
      function(read_length, asite_disp) {
        CalcAsiteFixedOneLength(
          reads_pos_length,
          min_read_length,
          read_length,
          asite_disp
        )
      }
    )
  if (colsum_out) {
    Asite_counts <- purrr::reduce(Asite_counts_bylength, `+`)
    return(Asite_counts)
  } else {
    # this has only as many columns as asite_disp_length,
    # probably LESS than data_mat
    Asite_counts_bylengthmat <- unlist(Asite_counts_bylength) %>%
      matrix(ncol = npos, byrow = TRUE)
    return(Asite_counts_bylengthmat)
  }
}

SumByFrame <- function(x, left, right) {
  # sum vector by 3nt frames 0,1,2
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  positions_frame0 <- seq(left, right, 3) # positions used to pick out frame 0 reads
  sums_byframe <- c(
    x[ positions_frame0 ] %>% sum(),
    x[ positions_frame0 + 1 ] %>% sum(),
    x[ positions_frame0 + 2 ] %>% sum()
  )
  return(sums_byframe)
}

SnapToCodon <- function(x, left, right, snapdisp=0L) {
  # snap nucleotide-aligned reads to codon position
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  #   snapdisp: integer any additional displacement in the snapping
  RcppRoll::roll_suml(x[(left:right) + snapdisp], n=3L, by=3L, fill = NULL)
}

NormByMean <- function(x,...) {
  x / mean(x,...)
}

GetGeneCodonPosReads1dsnap <- function(gene, dataset, hdf5file, left, right, 
                         min_read_length, 
                         asite_disp_length = data.frame(
                             read_length = c(28, 29, 30),
                             asite_disp = c(15, 15, 15)
                           ), 
                         snapdisp=0L) {
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hdf5file) # Get the matrix of read counts
  reads_asitepos <- CalcAsiteFixed(
    reads_pos_length, min_read_length,
    asite_disp_length
  )
  SnapToCodon(reads_asitepos,left,right,snapdisp)
}

GatherByFrameCodon <- function(x, left, right) {
  # gather vector by 3nt frames 0,1,2 for each position
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  positions_frame0 <- seq(left, right, 3) # positions to pick out frame 0 reads
  tibble(
    CodonPos = seq_len(length(positions_frame0)),
    Ct_fr0 = x[ positions_frame0 ],
    Ct_fr1 = x[ positions_frame0 + 1 ],
    Ct_fr2 = x[ positions_frame0 + 2 ]
  )
}

combinePValuesFisher <- function(p) {
  # Fisher's method (1-sided) to combine p-values
  pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = FALSE)
}

combinePValuesStouffer <- function(p) {
  # Stouffer's “inverse normal” method (1-sided) to combine p-values
  pnorm(sum(qnorm(p)) / sqrt(length(p)))
}

WilcoxTestFrame <- function(x, left, right) {
  # Wilcoxon rank-sum paired test that frame 0 has more reads
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  gathered_by_frame <- GatherByFrameCodon(x, left, right)

  wtresults_fr0vs1 <-
    wilcox.test(
      x = gathered_by_frame$Ct_fr0,
      y = gathered_by_frame$Ct_fr1,
      alternative = "greater", paired = TRUE, exact = FALSE
    )
  wtresults_fr0vs2 <-
    wilcox.test(
      x = gathered_by_frame$Ct_fr0,
      y = gathered_by_frame$Ct_fr2,
      alternative = "greater", paired = TRUE, exact = FALSE
    )

  return(c(
    pval_fr0vs1 = wtresults_fr0vs1$p.value,
    pval_fr0vs2 = wtresults_fr0vs2$p.value,
    pval_fr0vsboth =
      combinePValuesStouffer(c(
        wtresults_fr0vs1$p.value,
        wtresults_fr0vs2$p.value
      ))
  ))
}

GetGeneReadFrame <- function(gene, dataset, hdf5file, left, right, min_read_length,
                             asite_disp_length = data.frame(
                               read_length = c(28, 29, 30),
                               asite_disp = c(15, 15, 15)
                             )) {
  # example from vignette:
  #   GetGeneReadFrame(hdf5file, "YAL003W", dataset, 251, 871, min_read_length)
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hdf5file)
  reads_asitepos <- CalcAsiteFixed(
    reads_pos_length, min_read_length,
    asite_disp_length
  )
  sum_by_frame <- SumByFrame(reads_asitepos, left, right)
  wt_frame <- WilcoxTestFrame(reads_asitepos, left, right)
  tibble(
    gene = gene,
    Ct_fr0 = sum_by_frame[1],
    Ct_fr1 = sum_by_frame[2],
    Ct_fr2 = sum_by_frame[3],
    pval_fr0vs1 = wt_frame[1],
    pval_fr0vs2 = wt_frame[2],
    pval_fr0vsboth = wt_frame[3]
  )
}

CalcReadFrameProportion <- function(read_frame_df) {
  # calculate read frame proportions from data frame with read frame counts,
  stopifnot(all(c("Ct_fr0", "Ct_fr1", "Ct_fr2") %in% names(read_frame_df)))
  read_frame_df %>%
    mutate(
      Ct_all = Ct_fr0 + Ct_fr1 + Ct_fr2,
      p_fr0 = Ct_fr0 / Ct_all,
      p_fr1 = Ct_fr1 / Ct_all,
      p_fr2 = Ct_fr2 / Ct_all
    ) %>%
    return()
}

BoxplotReadFrameProportion <- function(read_frame_df, feat_names = "gene") {
  # Plot proportion of read frames as boxplot.
  rf_prop_long <- read_frame_df %>%
    CalcReadFrameProportion() %>%
    select(c(feat_names, "p_fr0", "p_fr1", "p_fr2")) %>%
    gather(-feat_names, key = "Frame", value = "Proportion") %>%
    mutate(Frame = factor(Frame,
      levels = c("p_fr0", "p_fr1", "p_fr2"),
      labels = 0:2
    ))
  ggplot(data = rf_prop_long, aes(x = Frame, colour = Frame, y = Proportion)) +
    geom_boxplot() +
    scale_y_continuous("Proportion, by feature", limits = c(0, 1), expand = c(0, 0)) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}


#####
## biases in nucleotide composition along mapped read lengths

GetNTReadPosition <- function(gene, dataset, hdf5file, lid, min_read_length) {
  reads_pos_len <- GetGeneDatamatrix(gene, dataset, hdf5file)[lid, ] # Get reads of a particular length
  reads_pos_len <- reads_pos_len[1:(length(reads_pos_len) - (lid + min_read_length - 1))] # Ignore reads whose 5' ends map close to the end of the 3' buffer
  pos <- rep(1:length(reads_pos_len), reads_pos_len) # nt positions weighted by number of reads mapping to it
  pos_IR <- IRanges::IRanges(start = pos, width = (lid + min_read_length - 1)) # Create an IRanges object for position-specific reads of a particular length
  return(pos_IR)
}

cons_mat <- function(gene, pos_IR, type = "count", cframe = 0, lid) {
  pos_IR_frame <- pos_IR[start(pos_IR) %% 3 == cframe] # Get position-specific reads of a particular length and ORF frame
  if (length(pos_IR_frame)) {
    pos_nt <- Biostrings::consensusMatrix(Biostrings::extractAt(coding_seqs[[gene]], pos_IR_frame))[1:4, ] # Get position-specific nucleotide counts
    if (type == "freq") {
      pos_nt <- pos_nt / colSums(pos_nt) # Select frequencies instead of counts
    }
  } else {
    pos_nt <- matrix(0, ncol = read_range[lid], nrow = 4, dimnames = list(c("A", "C", "G", "T"), NULL))
  }
  return(pos_nt)
}

comb_freq <- function(allfr) {
  alph <- c("A", "C", "G", "T")
  nt_sp_freq <- c()
  for (i in alph) {
    nt_sp_freq <- rbind(nt_sp_freq, colSums(allfr[rownames(allfr) == i, ])) # Get position-specific counts/freq of a nt across all genes for reads of a poarticular length and frame
  }
  nt_sp_freq <- t(nt_sp_freq / colSums(nt_sp_freq)) # Convert total counts to freq OR freq to normalized frequencies
  colnames(nt_sp_freq) <- alph
  return(nt_sp_freq)
}

#####
# read length and density functions

# read length-specific total counts stored as attributes of 'reads_total' in H5 file
GetGeneLength <- function(gene, dataset, hdf5file) {
  start_codon_pos <- rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "start_codon_pos"))[1]
  stop_codon_pos <- rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "stop_codon_pos"))[1]
  return(stop_codon_pos - start_codon_pos)
}

GetGeneReadsTotal <- function(gene, dataset, hdf5file) {
  rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "reads_total"))
}

GetGeneReadDensity <- function(gene, dataset, hdf5file, buffer = 50) {
  # buffer
  GetGeneReadsTotal(gene, dataset, hdf5file) / (GetGeneLength(gene, dataset, hdf5file) + 50)
}
