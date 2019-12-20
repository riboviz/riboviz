suppressMessages(library(getopt, quietly=T))
# Determine location of provenance.R relative to current file
source(file.path(dirname(getopt::get_Rscript_filename()), "provenance.R"))
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

# set ggplot2 theme for plots drawn after this; use dark on light theme
ggplot2::theme_set(theme_bw())

# define input options for optparse package
option_list <- list(
  make_option("--output-dir",
    type = "character", default = "./",
    help = "Output directory"
  ),
  make_option("--orf-fasta-file",
    type = "character", default = FALSE,
    help = "FASTA file with nt seq"
  ),
  make_option("--orf-gff-file",
    type = "character", default = NA,
    help = "riboviz generated GFF2/GFF3 annotation file"
  ),
  make_option("--num-processes",
    type = "integer", default = 1,
    help = "Number of cores for parallelization"
  ),
  make_option("--min-read-length",
    type = "integer", default = 10,
    help = "Minimum read length in H5 output"
  ),
  make_option("--max-read-length",
    type = "integer", default = 50,
    help = "Maximum read length in H5 output"
  ),
  make_option("--buffer",
    type = "integer", default = 250,
    help = "Length of flanking region around the CDS"
  ),
  make_option("--primary-id",
    type = "character", default = "gene_id",
    help = "Primary gene IDs to access the data (YAL001C, YAL003W, etc.)"
  ),
  make_option("--dataset",
    type = "character", default = "vignette",
    help = "Name of the dataset"
  ),
  make_option("--rpf",
    type = "logical", default = TRUE,
    help = "Is the dataset an RPF or mRNA dataset?"
  ),
  make_option("--features-file",
    type = "character", default = NA,
    help = "features file, columns are gene features and rows are genes"
  ),
  make_option("--do-pos-sp-nt-freq",
    type = "logical", default = TRUE,
    help = "do calculate the position-specific nucleotide frequency"
  ),
  make_option("--t-rna-file",
    type = "character", default = NA,
    help = "tRNA estimates in .tsv file"
  ),
  make_option("--codon-positions-file",
    type = "character", default = NA,
    help = "Codon positions in each gene in .Rdata file"
  ),
  make_option("--count-threshold",
    type = "integer", default = 64,
    help = "threshold for count of reads per gene to be included in plot"
  ),
  make_option("--output-prefix",
    type = "character", default = "",
    help = "Prefix for output files"
  ),
  make_option("--hd-file",
    type = "character", default = "output.h5",
    help = "Location of H5 output file"
  ),
  make_option("--nnt-buffer",
    type = "integer", default = 25,
    help = "n nucleotides of UTR buffer to include in metagene plots"
  ),
  make_option("--nnt-gene",
    type = "integer", default = 50,
    help = "n nucleotides of gene to include in metagene plots"
  ),
  make_option("--asite-disp-length-file",
    type = "character", default = NA,
    help = "asite displacement file
    table with one displacement per read length"
  )
)

print(get_version(get_Rscript_filename()))
# read in commandline arguments
opt <- optparse::parse_args(OptionParser(option_list = option_list),
                            convert_hyphens_to_underscores=TRUE)

attach(opt)

print("generate_stats_figs.R running with parameters:")
opt

# prepare files, opens hdf5 file connection
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

# read in positions of all exons/genes in GFF format and subset CDS locations
gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name

# read in coding sequences
coding_seqs <- readDNAStringSet(orf_fasta_file)

# range of read lengths between parameters set in config file
read_range <- min_read_length:max_read_length

# read in positions of all exons/genes in GFF format and subset CDS locations
gff <- readGFFAsGRanges(orf_gff_file)
gff_df <- gff %>% data.frame %>% as_tibble # @ewallace: tidy tibble version

# check for 3nt periodicity
print("Starting: Check for 3nt periodicity globally")

# function to get data matrix of read counts for gene and dataset from hdf5file
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

# function to get matrix of read counts from n_buffer before start codon to nnt_gene after
# for gene and dataset from hd5 file hdf5file, using UTR5 annotations in gff
GetGeneDatamatrix5start <- function(gene, dataset, hdf5file, gff,
                                    n_buffer = nnt_buffer,
                                    nnt_gene = nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  # @ewallace: replace this by gff_df?
  n_utr5 <- BiocGenerics::width(gff[gff$type == "UTR5" & gff$Name == gene])
  # if n_buffer bigger than length n_utr5, pad with zeros:
  if (n_utr5 >= n_buffer) {
    # if length n_utr5 bigger than n_buffer
    n_left5 <- n_utr5 - n_buffer + 1 # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  } else {
    # if length n_utr5 less than n_buffer
    n_left5 <- 1 # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = (n_buffer - n_utr5))
  }
  n_right3 <- n_utr5 + nnt_gene # column to end with (3'end)
  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  return(cbind(zeropad5_mat, data_mat_5start))
}

GetGeneDatamatrix3end <- function(gene, dataset, hdf5file, gff,
                                  n_buffer = nnt_buffer, nnt_gene = nnt_gene) {
  # get data matrix of read counts from nnt_gene before stop codon to n_buffer after
  # for gene and dataset from hd5 file hdf5file, using UTR3 annotations in gff
  # if n_buffer bigger than length n_utr3, pad with zeros.
  # CHECK startpos/off-by-one
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_all <- ncol(data_mat_all)
  # @ewallace: replace this by gff_df?
  n_utr3 <- BiocGenerics::width(gff[gff$type == "UTR3" & gff$Name == gene])
  n_left5 <- n_all - n_utr3 - nnt_gene + 1 # column to start from (5'end)
  if (n_utr3 >= n_buffer) {
    # length n_utr3 bigger than n_buffer
    n_right3 <- n_all - n_utr3 + n_buffer # column to end with (3'end)
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

GetNTPeriod <- function(hdf5file, gene, dataset, left, right) {
  # previous version of script; not currently used
  data_mat <- GetGeneDatamatrix(gene, dataset, hdf5file)
  pos_sum <- colSums(data_mat) # Position-specific sum of reads of all lengths
  pos_sum <- pos_sum[left:(ncol(data_mat) - right)] # Ignore reads in parts of the buffer region defined by left/right
  return(pos_sum)
}

# NOTE: Do not use mclapply when accessing H5 data

# get gene and position specific total counts for all read lengths
gene_poslen_counts_5start <-
  lapply(gene_names,
    GetGeneDatamatrix5start,
    dataset,
    hdf5file,
    gff,
    n_buffer = nnt_buffer,
    nnt_gene
  ) %>%
  Reduce("+", .) # Reduce is a BiocGenerics, there is a purrr alternative

# ribogrid_5start
gene_poslen_counts_5start %>%
  TidyDatamatrix(startpos = -nnt_buffer + 1, startlen = min_read_length) %>%
  plot_ribogrid() %>%
  ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogrid.pdf")),
    width = 6, height = 3
  )

gene_poslen_counts_5start %>%
  TidyDatamatrix(startpos = -nnt_buffer + 1, startlen = min_read_length) %>%
  barplot_ribogrid() %>%
  ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogridbar.pdf")),
    width = 6, height = 5
  )

gene_poslen_counts_3end <-
  lapply(gene_names,
    GetGeneDatamatrix3end,
    dataset,
    hdf5file,
    gff,
    n_buffer = nnt_buffer,
    nnt_gene
  ) %>%
  Reduce("+", .) # Reduce is a BiocGenerics, there is a purrr alternative

# summarize by adding different read lengths
gene_pos_counts_5start <- gene_poslen_counts_5start %>%
  TidyDatamatrix(startpos = -nnt_buffer + 1, startlen = min_read_length) %>%
  group_by(Pos) %>%
  summarize(Counts = sum(Counts))

gene_pos_counts_3end <- gene_poslen_counts_3end %>%
  TidyDatamatrix(startpos = -nnt_gene + 1, startlen = min_read_length) %>%
  group_by(Pos) %>%
  summarize(Counts = sum(Counts))

gene_pos_counts_bothends <- bind_rows(
  gene_pos_counts_5start %>% mutate(End = "5'"),
  gene_pos_counts_3end %>% mutate(End = "3'")
) %>%
  mutate(End = factor(End, levels = c("5'", "3'")))

# Plot
nt_period_plot <- ggplot(
  gene_pos_counts_bothends,
  aes(x = Pos, y = Counts)) +
  geom_line() +
  facet_wrap(~End, scales = "free") +
  labs(x = "Nucleotide Position", y = "Read counts")

# Save plot and file
ggsave(nt_period_plot, filename = file.path(output_dir, paste0(output_prefix, "3nt_periodicity.pdf")))
write.table(
  gene_pos_counts_bothends,
  file = file.path(output_dir, paste0(output_prefix, "3nt_periodicity.tsv")),
  sep = "\t",
  row = F,
  col = T,
  quote = F)
print("Completed: Check for 3nt periodicity globally")

## distribution of lengths of all mapped reads
print("Starting: Distribution of lengths of all mapped reads")

# read length-specific read counts stored as attributes of 'reads' in H5 file
gene_sp_read_length <- lapply(gene_names, function(x) {
  rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", x, "/", dataset, "/reads")), "reads_by_len"))
})


# sum reads of each length across all genes
read_length_data <- data.frame(
  Length = read_range,
  Counts = colSums(matrix(unlist(gene_sp_read_length), ncol = length(read_range), byrow = T))
)

# plot read lengths with counts
read_len_plot <- ggplot(read_length_data, aes(x = Length, y = Counts)) +
  geom_bar(stat = "identity")

# save read lengths plot and file
ggsave(read_len_plot, filename = file.path(output_dir, paste0(output_prefix, "read_lengths.pdf")))
write.table(
  read_length_data,
  file = file.path(output_dir, paste0(output_prefix, "read_lengths.tsv")),
  sep = "\t",
  row = F,
  col = T,
  quote = F
)

print("Completed: Distribution of lengths of all mapped reads")

## biases in nucleotide composition along mapped read lengths

print("Starting: Biases in nucleotide composition along mapped read lengths")

GetNTReadPosition <- function(hdf5file, gene, dataset, lid, min_read_length) {
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

if (do_pos_sp_nt_freq) {
  # This is in a conditional loop because it fails for some inputs
  # and has not been debugged.
  all_out <- c()
  for (lid in 1:length(read_range)) {
    out <- lapply(gene_names, function(x) {
      GetNTReadPosition(hdf5file, gene = as.character(x), dataset, lid = lid, min_read_length = min_read_length) # For each read length convert reads to IRanges
    })
    names(out) <- gene_names

    # Get position-specific nucleotide counts for reads in each frame
    fr0 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 0, lid = lid)
    }, mc.cores = num_processes)
    allfr0 <- do.call(rbind, fr0)
    fr1 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 1, lid = lid)
    }, mc.cores = num_processes)
    allfr1 <- do.call(rbind, fr1)
    fr2 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 2, lid = lid)
    }, mc.cores = num_processes)
    allfr2 <- do.call(rbind, fr2)

    # Get position-specific freq for all nts
    cnt_fr0 <- signif(comb_freq(allfr0), 3)
    cnt_fr1 <- signif(comb_freq(allfr1), 3)
    cnt_fr2 <- signif(comb_freq(allfr2), 3)

    output <- data.frame(rbind(cnt_fr0, cnt_fr1, cnt_fr2))
    all_out <- rbind(all_out, output)
  }

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

  all_out <- cbind(Length, Position, Frame, all_out)
  all_out[is.na(all_out)] <- 0

  # save file
  write.table(all_out, file = file.path(output_dir, paste0(output_prefix, "pos_sp_nt_freq.tsv")), sep = "\t", row = F, col = T, quote = F)

  print("Completed nucleotide composition bias table")
}

## calculate read frame for every annotated ORF

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

GatherByFrameCodon <- function(x, left, right) {
  # gather vector by 3nt frames 0,1,2 for each position
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  positions_frame0 <- seq(left, right, 3) # positions to pick out frame 0 reads
  tibble(
    CodonPos = 1:length(positions_frame0),
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

GetGeneReadFrame <- function(hdf5file, gene, dataset, left, right, min_read_length,
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

if (!is.na(asite_disp_length_file)) {
  print("Starting: Check for 3nt periodicity (frame) by Gene")
  asite_disp_length <- readr::read_tsv(asite_disp_length_file,
    comment = "#"
  )
  gene_read_frames <- gff_df %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::select(gene = seqnames, left = start, right = end) %>%
    purrr::pmap_dfr(GetGeneReadFrame,
      hdf5file = hdf5file,
      dataset = dataset,
      min_read_length = min_read_length,
      asite_disp_length = asite_disp_length
    )
  write.table(
    gene_read_frames,
    file = file.path(output_dir, paste0(output_prefix, "3ntframe_bygene.tsv")),
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )

  gene_read_frame_plot <- gene_read_frames %>%
    filter(Ct_fr0 + Ct_fr1 + Ct_fr2 > count_threshold) %>%
    BoxplotReadFrameProportion()

  # save read lengths plot and file
  ggsave(gene_read_frame_plot,
    filename = file.path(output_dir, paste0(output_prefix, "3ntframe_propbygene.pdf")),
    width = 3, height = 3
  )

  print("Completed: Check for 3nt periodicity (frame) by Gene")
}

## position specific distribution of reads

print("Starting: Position specific distribution of reads")

# codon-specific reads for RPF datasets
GetCodonPositionReads <- function(hdf5file, gene, dataset, left, right, min_read_length) {
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
GetMRNACoverage <- function(hdf5file, gene, dataset, left, right, read_range, min_read_length, buffer) {
  reads_pos <- GetGeneDatamatrix(gene, dataset, hdf5file) # Get the matrix of read counts
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

# For RPF datasets, generate codon-based position-specific reads
if (rpf) {
  # create empty matrix to store position-specific read counts
  out5p <- matrix(NA, nrow = length(gene_names), ncol = 500) # 5'
  out3p <- matrix(NA, nrow = length(gene_names), ncol = 500) # 3'

  out <- lapply(gene_names, function(gene) {
    GetCodonPositionReads(
      hdf5file,
      gene,
      dataset,
      left = (buffer - 15),
      right = (buffer + 11),
      min_read_length = min_read_length
      )
  }) # Get codon-based position-specific reads for each gene
  names(out) <- gene_names

  cc <- 1
  for (gene in gene_names) {
    tmp <- out[[gene]]
    # Only consider genes with at least count_threshold mapped reads along its CDS
    if (sum(tmp) >= count_threshold) {
      tmp <- tmp / mean(tmp)
      if (length(tmp) > 500) {
        out5p[cc, ] <- tmp[1:500]
        out3p[cc, ] <- rev(tmp)[1:500]
      } else {
        out5p[cc, 1:length(tmp)] <- tmp
        out3p[cc, 1:length(tmp)] <- rev(tmp)
      }
    }
    cc <- cc + 1
  }

  # Estimate position-specific mean and std error of mapped read counts
  m5p <- signif(apply(out5p, 2, mean, na.rm = T), 4)
  m3p <- signif(apply(out3p, 2, mean, na.rm = T), 4)
  s5p <- signif(apply(out5p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)
  s3p <- signif(apply(out3p, 2, function(x) {
    sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
  }), 4)

  # Normalize reads to last 50 codons of the 500-codon window.
  # This allows easy comparison between datasets
  s5p <- s5p / mean(m5p[450:500])
  s3p <- s3p / mean(m3p[450:500])
  m5p <- m5p / mean(m5p[450:500])
  m3p <- m3p / mean(m3p[450:500])

  # Create a dataframe to store the output for plots/analyses
  pos_sp_rpf_norm_reads <- data.frame(
    Position = c(1:500, 0:-499),
    Mean = c(m5p, m3p),
    SD = c(s5p, s3p),
    End = factor(rep(c("5'", "3'"), each = 500), levels = c("5'", "3'"))
  )

  # Plot
  pos_sp_rpf_norm_reads_plot <- ggplot(
    pos_sp_rpf_norm_reads,
    aes(Position, Mean, col = End)
  ) +
    geom_line() +
    facet_grid(~End, scales = "free") +
    guides(col = FALSE)

  # Save plot and file
  ggsave(pos_sp_rpf_norm_reads_plot, filename = file.path(output_dir, paste0(output_prefix, "pos_sp_rpf_norm_reads.pdf")))
  write.table(
    pos_sp_rpf_norm_reads,
    file = file.path(output_dir, paste0(output_prefix, "pos_sp_rpf_norm_reads.tsv")),
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
}

# for mRNA datasets, generate nt-based position-specific reads
if (!rpf) {
  # create empty matrix to store position-specific read counts
  out5p <- matrix(NA, nrow = length(gene_names), ncol = 1500) # 5'
  out3p <- matrix(NA, nrow = length(gene_names), ncol = 1500) # 3'

  out <- lapply(gene_names, function(gene) {
    GetMRNACoverage(
      hdf5file,
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

  # plot
  pos_sp_mrna_norm_coverage_plot <- ggplot(pos_sp_mrna_norm_coverage, aes(Position, Mean, col = End)) +
    geom_line() +
    facet_grid(~End, scales = "free") +
    guides(col = FALSE)

  # Save plot and file
  ggsave(pos_sp_mrna_norm_coverage_plot, filename = file.path(output_dir, paste0(output_prefix, "pos_sp_mrna_norm_coverage.pdf")))
  write.table(
    pos_sp_mrna_norm_coverage,
    file = file.path(output_dir, paste0(output_prefix, "pos_sp_mrna_norm_coverage.tsv")),
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
}

print("Completed: Position specific distribution of reads")

## Calculate TPMs of genes

print("Starting: Calculate TPMs of genes")

# read length-specific total counts stored as attributes of 'reads_total' in H5 file
GetGeneLength <- function(gene, hdf5file) {
  start_codon_pos <- rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "start_codon_pos"))[1]
  stop_codon_pos <- rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "stop_codon_pos"))[1]
  return(stop_codon_pos - start_codon_pos)
}

GetGeneReadsTotal <- function(gene, hdf5file) {
  rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "reads_total"))
}

GetGeneReadDensity <- function(gene, hdf5file, buffer = 50) {
  # buffer
  GetGeneReadsTotal(gene, hdf5file) / (GetGeneLength(gene, hdf5file) + 50)
}

# calculate transcripts per million (TPM)
gene_sp_reads <- sapply(gene_names, GetGeneReadsTotal, hdf5file)
reads_per_b <- sapply(gene_names, GetGeneReadDensity, hdf5file)

tpms <- data.frame(
  ORF = gene_names,
  readcount = gene_sp_reads,
  rpb = reads_per_b,
  tpm = reads_per_b * 1e6 / sum(reads_per_b)
)

# write out to *_tpms.tsv
write.table(
  tpms,
  file = file.path(output_dir, paste0(output_prefix, "tpms.tsv")),
  sep = "\t",
  row = F,
  col = T,
  quote = F
)


## Correlations between TPMs of genes with their sequence-based features

# Correlate TPMs of genes with sequence-based features, skip if missing features_file
if (!is.na(features_file)) {
  print("Starting: Correlations between TPMs of genes with their sequence-based features")

  features <- read.table(features_file, h = T)

  # Prepare data for plot
  # Consider only genes with at least count_threshold mapped reads
  features_plot_data <- merge(features, tpms, by = "ORF") %>%
    filter(readcount >= count_threshold, !is.na(ORF)) %>%
    select(-readcount, -rpb) %>%
    gather(Feature, Value, -ORF, -tpm)

  features_plot <- ggplot(features_plot_data, aes(x = tpm, y = Value)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~Feature, scales = "free") +
    scale_x_log10() +
    geom_smooth(method = "lm") +
    xlab("TPM (transcripts per million)")

  # Save plot and file
  ggsave(features_plot, filename = file.path(output_dir, paste0(output_prefix, "features.pdf")))

  print("Completed: Correlations between TPMs of genes with their sequence-based features")
} else {
  print("Skipped: Correlations between TPMs of genes with their sequence-based features - features_file.tsv not provided")
}

## Codon-specific ribosome densities for correlations with tRNAs

# Codon-specific ribosome density for tRNA correlation; skip if missing t_rna_file & codon_positions_file
if (!is.na(t_rna_file) & !is.na(codon_positions_file)) {
  print("Starting: Codon-specific ribosome densities for correlations with tRNAs")

  # Only for RPF datasets
  if (rpf) {
    # This still depends on yeast-specific arguments and should be edited.
    yeast_tRNAs <- read.table(t_rna_file, h = T) # Read in yeast tRNA estimates
    load(codon_positions_file) # Position of codons in each gene (numbering ignores first 200 codons)
    # Reads in an object named "codon_pos"
    out <- lapply(gene_names, function(gene) {
      # From "Position specific distribution of reads" plot
      GetCodonPositionReads(hdf5file, gene, dataset, left = (buffer - 15), right = (buffer + 11), min_read_length = min_read_length)
    }) # Get codon-based position-specific reads for each gene
    names(out) <- gene_names

    gene_len <- sapply(out, length) # Calculate gene length in codons
    out <- out[gene_len > 201] # Ignore genes with <=200 sense codons

    trim_out <- lapply(out, function(x) {
      x[201:(length(x) - 1)]
    }) # Trim first 200 codons and stop codon from each gene
    read_counts_trim <- sapply(trim_out, sum) # Calculate read counts in trimmed genes
    trim_out <- trim_out[read_counts_trim >= count_threshold] # Ignore genes with fewer than count_threshold mapped reads

    norm_out <- lapply(trim_out, function(x) {
      x / mean(x)
    }) # Normalize reads in each gene by their mean

    # Calculate codon-specific mean ribosome-densities at A/P/E sites of the mapped reads
    a_mn <- sapply(names(codon_pos), function(codon) {
      mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
        pos <- as.numeric(a[2])
        norm_out[[a[1]]][pos]
      })), na.rm = T)
    })
    p_mn <- sapply(names(codon_pos), function(codon) {
      mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
        pos <- as.numeric(a[2]) + 1
        norm_out[[a[1]]][pos]
      })), na.rm = T)
    })
    e_mn <- sapply(names(codon_pos), function(codon) {
      mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
        pos <- as.numeric(a[2]) + 2
        norm_out[[a[1]]][pos]
      })), na.rm = T)
    })

    # Sort the values
    A <- a_mn[order(names(codon_pos))]
    P <- p_mn[order(names(codon_pos))]
    E <- e_mn[order(names(codon_pos))]

    cod_dens_tRNA <- cbind(yeast_tRNAs, A, P, E)

    # Prepare data for plot
    cod_dens_tRNA_wide <- cod_dens_tRNA %>%
      gather(tRNA_type, tRNA_value, 3:6) %>%
      gather(Site, Ribodens, 3:5)

    cod_dens_tRNA_plot <- ggplot(cod_dens_tRNA_wide, aes(x = tRNA_value, y = Ribodens)) +
      geom_point(alpha = 0.3) +
      facet_grid(Site ~ tRNA_type, scales = "free_x") +
      geom_smooth(method = "lm") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Save plot and file
    ggsave(cod_dens_tRNA_plot, filename = file.path(output_dir, paste0(output_prefix, "codon_ribodens.pdf")))
    write.table(
      cod_dens_tRNA,
      file = file.path(output_dir, paste0(output_prefix, "codon_ribodens.tsv")),
      sep = "\t",
      row = F,
      col = T,
      quote = F
    )
  }

  print("Completed: Codon-specific ribosome densities for correlations with tRNAs")
} else {
  print("Skipped: Codon-specific ribosome densities for correlations with tRNAs - t-rna-file and/or codon-positions-file not provided")
}
