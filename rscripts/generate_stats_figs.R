# generate_stats_figs.R
# currently (17 December 2019) undergoing refactoring by EW
# TODO:
# - tidy get_cod_pos
# - functionalize gene_poslen_counts, wrapper functions
# - tidy gene_sp_read_length to GetGeneReadLength
# - tidy do_pos_sp_nt_freq
# - functionalize gene_read_frames
# - tidy GetCodonPositionReads wrapper
# - functionalize pos_sp_rpf_norm_reads
# - tidy GetMRNACoverage wrapper
# - functionalize pos_sp_mrna_norm_coverage
# - functionalize GetTPMs into *_tpms.tsv
# - refactor Codon-specific ribosome densities
# - replace write.table with write_tsv
# - duplicate generate_stats_figs.R as .Rmd file
# - replace Reduce, lapply, sapply by purrr functions?
#
## Should use relative path for rscripts??
source("rscripts/read_count_functions.R")

# set ggplot2 theme for plots drawn after this; use dark on light theme
ggplot2::theme_set(theme_bw())

# define input options for optparse package
option_list <- list(
  make_option("--dir_out",
    type = "character", default = "./",
    help = "Output directory"
  ),
  make_option("--orf_fasta",
    type = "character", default = FALSE,
    help = "FASTA file with nt seq"
  ),
  make_option("--orf_gff_file",
    type = "character", default = NA,
    help = "riboviz generated GFF2/GFF3 annotation file"
  ),
  make_option("--Ncores",
    type = "integer", default = 1,
    help = "Number of cores for parallelization"
  ),
  make_option("--MinReadLen",
    type = "integer", default = 10,
    help = "Minimum read length in H5 output"
  ),
  make_option("--MaxReadLen",
    type = "integer", default = 50,
    help = "Maximum read length in H5 output"
  ),
  make_option("--Buffer",
    type = "integer", default = 250,
    help = "Length of flanking region around the CDS"
  ),
  make_option("--PrimaryID",
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
  make_option("--features_file",
    type = "character", default = NA,
    help = "features file, columns are gene features and rows are genes"
  ),
  make_option("--do_pos_sp_nt_freq",
    type = "logical", default = TRUE,
    help = "do calculate the position-specific nucleotide frequency"
  ),
  make_option("--t_rna",
    type = "character", default = NA,
    help = "tRNA estimates in .tsv file"
  ),
  make_option("--codon_pos",
    type = "character", default = NA,
    help = "Codon positions in each gene in .Rdata file"
  ),
  make_option("--count_threshold",
    type = "integer", default = 64,
    help = "threshold for count of reads per gene to be included in plot"
  ),
  make_option("--out_prefix",
    type = "character", default = "out",
    help = "Prefix for output files"
  ),
  make_option("--hdFile",
    type = "character", default = "output.h5",
    help = "Location of H5 output file"
  ),
  make_option("--nnt_buffer",
    type = "integer", default = 25,
    help = "n nucleotides of UTR buffer to include in metagene plots"
  ),
  make_option("--nnt_gene",
    type = "integer", default = 50,
    help = "n nucleotides of gene to include in metagene plots"
  ),
  make_option("--asite_disp_length_file",
    type = "character", default = NA,
    help = "asite displacement file
    table with one displacement per read length"
  )
)

#####
# read in commandline arguments
opt <- optparse::parse_args(OptionParser(option_list = option_list))
attach(opt)

print("generate_stats_figs.R running with parameters:")
opt

# prepare files, opens hdf5 file connection
hdf5file <- rhdf5::H5Fopen(hdFile) # filehandle for the h5 file

# read in positions of all exons/genes in GFF format and subset CDS locations
gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name

# read in coding sequences
coding_seqs <- readDNAStringSet(orf_fasta)

# range of read lengths between parameters set in config file
read_range <- MinReadLen:MaxReadLen

# read in positions of all exons/genes in GFF format and subset CDS locations
gff <- readGFFAsGRanges(orf_gff_file)
gff_df <- gff %>% data.frame %>% as_tibble # @ewallace: tidy tibble version

#####
# check for 3nt periodicity
print("Starting: Check for 3nt periodicity globally")

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
  TidyDatamatrix(startpos = -nnt_buffer + 1, startlen = MinReadLen) %>%
  plot_ribogrid() %>%
  ggsave(
    filename = paste0(out_prefix, "_startcodon_ribogrid.pdf"),
    width = 6, height = 3
  )

gene_poslen_counts_5start %>%
  TidyDatamatrix(startpos = -nnt_buffer + 1, startlen = MinReadLen) %>%
  barplot_ribogrid() %>%
  ggsave(
    filename = paste0(out_prefix, "_startcodon_ribogridbar.pdf"),
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
  TidyDatamatrix(startpos = -nnt_buffer + 1, startlen = MinReadLen) %>%
  group_by(Pos) %>%
  summarize(Counts = sum(Counts))

gene_pos_counts_3end <- gene_poslen_counts_3end %>%
  TidyDatamatrix(startpos = -nnt_gene + 1, startlen = MinReadLen) %>%
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
ggsave(nt_period_plot, filename = paste0(out_prefix, "_3nt_periodicity.pdf"))
write.table(
  gene_pos_counts_bothends,
  file = paste0(out_prefix, "_3nt_periodicity.tsv"),
  sep = "\t",
  row = F,
  col = T,
  quote = F)
print("Completed: Check for 3nt periodicity globally")

#####
## distribution of lengths of all mapped reads
print("Starting: Distribution of lengths of all mapped reads")

# read length-specific read counts stored as attributes of 'reads' in H5 file
gene_sp_read_length <- lapply(gene_names, function(gene) {
  # TODO: GetGeneReadLength(gene,dataset,hdf5file)
  rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "reads_by_len"))
})


# sum reads of each length across all genes
read_length_data <- data.frame(
  Length = read_range,
  # make this tidy
  Counts = colSums(matrix(unlist(gene_sp_read_length), 
                          ncol = length(read_range), byrow = T))
)

# plot read lengths with counts
read_len_plot <- ggplot(read_length_data, aes(x = Length, y = Counts)) +
  geom_bar(stat = "identity")

# save read lengths plot and file
ggsave(read_len_plot, filename = paste0(out_prefix, "_read_lengths.pdf"))
write.table(
  read_length_data,
  file = paste0(out_prefix, "_read_lengths.tsv"),
  sep = "\t",
  row = F,
  col = T,
  quote = F
)

print("Completed: Distribution of lengths of all mapped reads")

#####
## biases in nucleotide composition along mapped read lengths

print("Starting: Biases in nucleotide composition along mapped read lengths")

if (do_pos_sp_nt_freq) {
  # This is in a conditional loop because it fails for some inputs
  # and has not been debugged. Needs to be rewritten in tidyverse
  all_out <- c()
  for (lid in seq_len(length(read_range))) {
    out <- lapply(gene_names, function(x) {
      # For each read length convert reads to IRanges
      GetNTReadPosition(gene = as.character(x), 
                        dataset=dataset,
                        hdf5file = hdf5file, 
                        lid = lid, MinReadLen = MinReadLen)
    })
    names(out) <- gene_names

    # Get position-specific nucleotide counts for reads in each frame
    fr0 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 0, lid = lid)
    }, mc.cores = Ncores)
    allfr0 <- do.call(rbind, fr0)
    fr1 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 1, lid = lid)
    }, mc.cores = Ncores)
    allfr1 <- do.call(rbind, fr1)
    fr2 <- mclapply(gene_names, function(gene) {
      cons_mat(gene = gene, pos_IR = out[[gene]], cframe = 2, lid = lid)
    }, mc.cores = Ncores)
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
  write.table(all_out, file = paste0(out_prefix, "_pos_sp_nt_freq.tsv"), sep = "\t", row = F, col = T, quote = F)

  print("Completed nucleotide composition bias table")
}

## calculate read frame for every annotated ORF


if (!is.na(asite_disp_length_file)) {
  print("Starting: Check for 3nt periodicity (frame) by Gene")
  asite_disp_length <- readr::read_tsv(asite_disp_length_file,
    comment = "#"
  )
  # TODO
  gene_read_frames <- gff_df %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::select(gene = seqnames, left = start, right = end) %>%
    purrr::pmap_dfr(GetGeneReadFrame,
      hdf5file = hdf5file,
      dataset = dataset,
      MinReadLen = MinReadLen,
      asite_disp_length = asite_disp_length
    )
  write.table(
    gene_read_frames,
    file = paste0(out_prefix, "_3ntframe_bygene.tsv"),
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
    filename = paste0(out_prefix, "_3ntframe_propbygene.pdf"),
    width = 3, height = 3
  )

  print("Completed: Check for 3nt periodicity (frame) by Gene")
}

#####
## position specific distribution of reads

print("Starting: Position specific distribution of reads")

# For RPF datasets, generate codon-based position-specific reads
if (rpf) {
  # create empty matrix to store position-specific read counts
  out5p <- matrix(NA, nrow = length(gene_names), ncol = 500) # 5'
  out3p <- matrix(NA, nrow = length(gene_names), ncol = 500) # 3'
  
  # TODO
  out <- lapply(gene_names, function(gene) {
    GetCodonPositionReads(
      gene,
      dataset,
      hdf5file,
      left = (Buffer - 15),
      right = (Buffer + 11),
      MinReadLen = MinReadLen
      )
  }) # Get codon-based position-specific reads for each gene
  names(out) <- gene_names

  cc <- 1
  for (gene in gene_names) {
    tmp <- out[[gene]]
    # Only consider genes with at least count_threshold mapped reads along its CDS
    # Here we should have a function taking the count_threshold and the length (500)
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
  # TODO: 50 should be a parameter
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
  ggsave(pos_sp_rpf_norm_reads_plot, filename = paste0(out_prefix, "_pos_sp_rpf_norm_reads.pdf"))
  write.table(
    pos_sp_rpf_norm_reads,
    file = paste0(out_prefix, "_pos_sp_rpf_norm_reads.tsv"),
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
  
  # TODO
  out <- lapply(gene_names, function(gene) {
    GetMRNACoverage(
      gene,
      dataset,
      hdf5file,
      left = (Buffer - 49),
      right = (Buffer - 3),
      MinReadLen = MinReadLen,
      read_range,
      Buffer
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
  ggsave(pos_sp_mrna_norm_coverage_plot, filename = paste0(out_prefix, "_pos_sp_mrna_norm_coverage.pdf"))
  write.table(
    pos_sp_mrna_norm_coverage,
    file = paste0(out_prefix, "_pos_sp_mrna_norm_coverage.tsv"),
    sep = "\t",
    row = F,
    col = T,
    quote = F
  )
}

print("Completed: Position specific distribution of reads")

## Calculate TPMs of genes

print("Starting: Calculate TPMs of genes")

# calculate transcripts per million (TPM)
# TODO: make this a single function
gene_sp_reads <- sapply(gene_names, GetGeneReadsTotal, dataset, hdf5file)
reads_per_b <- sapply(gene_names, GetGeneReadDensity, dataset, hdf5file)

tpms <- data.frame(
  ORF = gene_names,
  readcount = gene_sp_reads,
  rpb = reads_per_b,
  tpm = reads_per_b * 1e6 / sum(reads_per_b)
)

# write out to *_tpms.tsv
write.table(
  tpms,
  file = paste0(out_prefix, "_tpms.tsv"),
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
  ggsave(features_plot, filename = paste0(out_prefix, "_features.pdf"))

  print("Completed: Correlations between TPMs of genes with their sequence-based features")
} else {
  print("Skipped: Correlations between TPMs of genes with their sequence-based features - features_file.tsv not provided")
}

## Codon-specific ribosome densities for correlations with tRNAs

# Codon-specific ribosome density for tRNA correlation; skip if missing t_rna & codon_pos
if (!is.na(t_rna) & !is.na(codon_pos)) {
  print("Starting: Codon-specific ribosome densities for correlations with tRNAs")

  # Only for RPF datasets
  if (rpf) {
    # TODO: This still depends on yeast-specific arguments and should be edited.
    yeast_tRNAs <- read.table(t_rna, h = T) # Read in yeast tRNA estimates
    load(codon_pos) # Position of codons in each gene (numbering ignores first 200 codons)
    # Reads in an object named "codon_pos"
    out <- lapply(gene_names, function(gene) {
      # From "Position specific distribution of reads" plot
      GetCodonPositionReads(gene, dataset, hdf5file, left = (Buffer - 15), right = (Buffer + 11), MinReadLen = MinReadLen)
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
    
    
    # TODO: figure this out
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
    ggsave(cod_dens_tRNA_plot, filename = paste0(out_prefix, "_codon_ribodens.pdf"))
    write.table(
      cod_dens_tRNA,
      file = paste0(out_prefix, "_codon_ribodens.tsv"),
      sep = "\t",
      row = F,
      col = T,
      quote = F
    )
  }

  print("Completed: Codon-specific ribosome densities for correlations with tRNAs")
} else {
  print("Skipped: Codon-specific ribosome densities for correlations with tRNAs - t_rna file and/or codon_pos file not provided")
}
