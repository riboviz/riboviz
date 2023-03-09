#' Convert BAM files to riboviz HDF5 files.
#'
#' Given a GFF file and a BAM file, this script creates some HDF5
#' files with information about a feature (e.g. CDS, ORF, or uORF).
#'
#' See `BamToH5` below.
#'
#' @export

suppressMessages(library(GenomeInfoDb, quietly = T))
suppressMessages(library(IRanges, quietly = T))
suppressMessages(library(Rsamtools, quietly = T))
suppressMessages(library(rtracklayer, quietly = T))
suppressMessages(library(rhdf5, quietly = T))
suppressMessages(library(parallel, quietly = T))

#' Compute matrix of read counts with start position and read length
#' on a single gene.
#'
#' @param gene_location Gene coordinates from GFF (GRanges).
#' @param bam_file BAM file (character).
#' @param read_lengths List of read lengths (integer).
#' @param left_flank Length of left-hand flanking region (integer).
#' @param right_flank Length of right-hand flanking region (integer).
#' @param mult_exon If `TRUE` use only the first supplied exon for
#' gene coordinates (logical).
#' @return matrix of integer counts for each start position (rows)
#' and read length (columns) (double).
#'
#' @export
ReadsToCountMatrix <- function(gene_location, bam_file, read_lengths,
  left_flank, right_flank, mult_exon = TRUE) {

  if (!mult_exon) {
    # If the gene is a single exon gene but multiple GRanges specified
    # use the first GRange.
    gene_location <- gene_location[1]
  }

  num_read_lengths <- length(read_lengths)
  gene_length <- sum(sum(IRanges::coverage(gene_location)))

  if (as.character(strand(gene_location)[1]) == "-") {
    # feature is on negative strand.
    tmp <- left_flank
    left_flank <- right_flank
    right_flank <- tmp
  }

  # Expand genomic locations to include flanking region positions.
  if (length(gene_location) == 1) {
    # Single exon gene.
    rtracklayer::start(gene_location) <-
      rtracklayer::start(gene_location) - left_flank
    rtracklayer::end(gene_location) <-
      rtracklayer::end(gene_location) + right_flank
  } else {
    # Multiple exon gene.
    rtracklayer::start(gene_location)[rtracklayer::start(gene_location$start) == min(rtracklayer::start(gene_location$start))] <-
      rtracklayer::start(gene_location$start)[rtracklayer::start(gene_location$start) == min(rtracklayer::start(gene_location$start))] - left_flank
    rtracklayer::end(gene_location)[rtracklayer::end(gene_location$end) == max(rtracklayer::end(gene_location$end))] <-
      rtracklayer::end(gene_location$end)[rtracklayer::end(gene_location$end) == max(rtracklayer::end(gene_location$end))] + right_flank
  }

  # Read gene's strand, pos, qwidth data from BAM file.
  bam_what <- c("strand", "pos", "qwidth")
  bam_param <- Rsamtools::ScanBamParam(
    which = gene_location, what = bam_what)
  bam_data <- Rsamtools::scanBam(bam_file, param = bam_param)

  # Subset reads that are on the same strand at the genomic location.
  read_strand <- unlist(lapply(bam_data, function(x) x$strand))
  read_location <- unlist(lapply(bam_data, function(x) x$pos))[
    read_strand == as.factor(strand(gene_location)[1])]
  read_width <- unlist(lapply(bam_data, function(x) x$qwid))[
    read_strand == as.factor(strand(gene_location)[1])]

  # Column numbers based on genomic position.
  nucleotide_pos <- unlist(which(coverage(gene_location)[
    GenomeInfoDb::seqnames(gene_location)[1]] == 1))

  # If the specified flanking regions of a gene end up outside the
  # chromosome locations (<0), add pseudo columns with negative
  # numbers.
  if (min(rtracklayer::start(gene_location)) < min(nucleotide_pos)) {
    nucleotide_pos <- c(min(rtracklayer::start(gene_location)):0,
                            nucleotide_pos)
  }

  # Count reads whose 5' ends map to each nucleotide.
  read_counts <- matrix(0,
                        nrow = num_read_lengths,
                        ncol = (gene_length + left_flank + right_flank))
  # Cast to be matrix of integer.
  mode(read_counts) <- "integer"
  if (all(strand(gene_location) == "+")) {
    # Genes on positive strands.
    row <- 1
    for (i in read_lengths) {
      # Subset reads of length i.
      read_loc_len_i <- read_location[read_width == i]
      # Count reads of length i.
      count_reads_len_i <- table(factor(read_loc_len_i,
                                        levels = nucleotide_pos))
      read_counts[row, ] <- c(count_reads_len_i)
      row <- row + 1
    }
  }
  if (all(strand(gene_location) == "-")) {
    # Genes on negative strands.
    row <- 1
    for (i in read_lengths) {
      # Subset reads of length i.
      read_loc_len_i <- read_width[read_width == i] + read_location[read_width == i] - 1
      # Count reads of length i.
      count_reads_len_i <- table(factor(read_loc_len_i,
                                        levels = nucleotide_pos))
      read_counts[row, ] <- c(rev(count_reads_len_i))
      row <- row + 1
    }
  }
  return(read_counts);
}

#' Convert BAM file to riboviz HDF5 file(s).
#'
#' @param bam_file BAM input file (character).
#' @param orf_gff_file GFF2/GFF3 Matched genome feature file,
#' specifying coding sequences locations (start and stop coordinates)
#' within the transcripts (GTF/GFF3 file) (character).
#' @param feature Feature e.g. `CDS`, `ORF`, or `uORF` (character).
#' @param min_read_length Minimum read length in H5 output (integer).
#' @param max_read_length Maximum read length in H5 output (integer).
#' @param buffer Length of flanking region around the feature. Used
#' only if `is_riboviz_gff` is `FALSE`. (logical). (integer).
#' @param primary_id Primary gene IDs to access the data (character).
#' @param secondary_id Secondary gene IDs to access the data (character).
#' @param dataset Human-readable name of the dataset (character).
#' @param stop_in_feature Are stop codons part of the feature
#' annotations in `orf_gff_file`? Used only if `is_riboviz_gff`
#' is `FALSE`. (logical).
#' @param is_riboviz_gff Does `orf_gff_file` contain 3 elements per
#' gene - UTR5, feature, and UTR3? (logical).
#' @param hd_file H5 output file (character).
#' @param num_processes Number of processes to parallelize over
#' (integer).
#'
#' All reads are mapped to their 5' ends.
#'
#' `primary_id` is the name of an attribute in `orf_gff_file`
#' (e.g. `Name`) expected to hold gene names.
#'
#' `secondary_id` is the name of an attribute in `orf_gff_file`
#' (e.g. `ID`) expected to hold gene names or, if none, `NULL`. If
#' provided then these alternative gene names are used to create
#' symbolic links in the H5 files to the entries for each gene.
#'
#' If `is_riboviz_gff` then:
#'
#' * Feature (e.g. `CDS`, `ORF`, or `uORF`), `UTR5` and `UTR3` entries
#'   from `orf_gff_file` are used.
#' * `buffer` is ignored.
#'
#' If not `is_riboviz_gff` then:
#'
#' * Feature (e.g. `CDS`, `ORF`, or `uORF`) entries from
#'   `orf_gff_file` are used.
#' * `UTR5` and `UTR3` entries from `orf_gff_file` are ignored.
#' * `buffer` is used as the width of left and right flanks.
#' * `stop_in_feature` states where the stop codon is located (true if
#'   within the feature, false otherwise).
#'
#' An HDF5 file, `hd_file`, is created that has external links to
#' complementary HDF5 files, named `<hd_file>.1`,  `<hd_file>.2`
#' etc. each of which hold the data for a subset of genes. The number
#' of data files depends on the number of processes
#' (`num_processes`). For example, if `num_processes` is 1 then the
#' output files will be `<hd_file>` (external links file) and
#' `<hd_file>.1` (data file). If `num_processes` is 4 then the
#' output files will be `<hd_file>` (external links file) and
#' `<hd_file>.1`, `<hd_file>.2`, `<hd_file>.3`, `<hd_file>.4`
#' (data files).  
#'
#' Each complementary HDF5 data file has approximately the same number
#' of genes i.e. approximately number of genes / `num_processes`.
#'
#' On output, `hd_file` has the following structure:
#'
#' ```
#' GROUP "/" {
#'   EXTERNAL_LINK "<gene<1>>" {
#'      TARGETFILE "<hd_file>.1"
#'      TARGETPATH "<gene<1>>"
#'   }
#'  ...
#'   EXTERNAL_LINK "<gene<m>>" {
#'      TARGETFILE "<hd_file>.2"
#'      TARGETPATH "<gene>m>>"
#'   }
#'  ...
#'   EXTERNAL_LINK "<gene<n>>" {
#'      TARGETFILE "<hd_file>.<num_processes>"
#'      TARGETPATH "<gene<n>>"
#'   }
#' ...
#' }
#' ```
#'
#' The complementary HDF5 data files, `<hd_file>.<i>`, each have the
#' following structure.
#'
#' The `reads` group, `/<gene>/<dataset>/reads`, for a `<gene>` has
#' several attributes associated with it. These are summary statistics
#' and other information about the gene and dataset within the `reads`
#' group. The list of attributes are as follows.
#'
#' | Attribute | Description | Origin |
#' | --------- |------------ | ------ |
#' | `buffer_left` | Number of nucleotides upstream of the start codon (ATG) (UTR5 length) | EITHER position of start codon (from `orf_gff_file`) - 1 OR, if `is_riboviz_gff` is false, `buffer` |
#' | `buffer_right` | Number of nucleotides downstream of the stop codon (TAA/TAG/TGA) (UTR3 length) | From rom `orf_gff_file` OR, if `is_riboviz_gff` is false, `buffer` + feature length from `orf_gff_file` + `buffer` - `stop_codon_pos[3]` (see below) |
#' | `start_codon_pos` | Positions corresponding to start codon of feature (typically 251, 252, 253) | From `orf_gff_file` |
#' | `stop_codon_pos` | Positions corresponding to stop codon of feature | From `orf_gff_file`. If `is_riboviz_gff` is false and `stop_in_feature` is false, position is last nucleotide of feature + 1 |
#' | `lengths` | Lengths of mapped reads | Range from `min_read_length` to `max_read_length` (e.g. 10,...,50) |
#' | `reads_by_len` | Counts of number of ribosome sequences of each length | BAM file |
#' | `reads_total` | Total number of ribosome sequences | BAM file. Equal to number of non-zero reads in `reads_by_len` |
#'
#' Each complementary HDF5 data file, `<hd_file>.<i>`, is organized in
#' a hierarchy, `/<gene>/<dataset>/reads/data`, where `<gene>` is each
#' gene name in `bam_file` and `orf_gff_file`.
#'
#' The `data` table in `/<gene>/<dataset>/reads/data`, for a `<gene>`
#' has the positions and lengths of ribosome sequences within the
#' organism data (determined from `bam_file`). It is an integer table
#' with each row representing a read length and columns representing
#' nucleotide positions. The first row corresponds to reads of length
#' `min_read_length` and the last row corresponds to reads of length
#' `max_read_length`.
#'
#' A template HDF5 data file, showing how the complementary HDF5 data
#' files relate to the arguments to this function, and the contents of
#' `orf_gff_file` and `bam_file` is as follows:
#'
#' ```
#' HDF5 "<hd_file>.<i>" {
#' GROUP "/" {
#'   GROUP "<gene>" {
#'     GROUP "<dataset>" {
#'       GROUP "reads" {
#'         ATTRIBUTE "buffer_left" {
#'            DATATYPE  H5T_STD_I32LE
#'            DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
#'            DATA {
#'              (0): <Feature UTR5 length from orf_gff_file OR buffer (if is_riboviz_gff is false)>
#'            }
#'         }
#'         ATTRIBUTE "buffer_right" {
#'            DATATYPE  H5T_STD_I32LE
#'            DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
#'            DATA {
#'              (0): <Feature UTR3 length from orf_gff_file OR buffer (if is_riboviz_gff is false)>
#'            }
#'         }
#'         ATTRIBUTE "lengths" {
#'            DATATYPE  H5T_STD_I32LE
#'            DATASPACE  SIMPLE { (<read_length>) / (<read_length>) }
#'            DATA {
#'              (0): <min_read_length>, <min_read_length> + 1, ... , <min_read_length> + <m> _ 1,
#'              (<m>): <min_read_length>+<m>, <min_read_length> + <m> + 1, ... , <min_read_length> + <n> _ 1,
#'              (<n>): <min_read_length>+<n>, <min_read_length> + <n> + 1, ... , <min_read_length> + <max_read_length>
#'            }
#'         }
#'         ATTRIBUTE "reads_by_len" {
#'           DATATYPE  H5T_STD_I32LE
#'           DATASPACE  SIMPLE { (<read_length>) / (<read_length>) }
#'           DATA {
#'             (0): <see below>
#'             (<m>): <see below>
#'             (<n>): <see below>
#'           }
#'         }
#'         ATTRIBUTE "reads_total" {
#'           DATATYPE  H5T_STD_I32LE
#'           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
#'           DATA {
#'             (0): <number of non-zero values in reads_by_len>
#'           }
#'         }
#'         ATTRIBUTE "start_codon_pos" {
#'           DATATYPE  H5T_STD_I32LE
#'           DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
#'           DATA {
#'             (0): <position of 1st nt of feature start codon from orf_gff_file>,
#'                  <position of 2nd nt of feature start codon from orf_gff_file>,
#'                  <position of 3rd nt of feature start codon from orf_gff_file>
#'           }
#'         }
#'         ATTRIBUTE "stop_codon_pos" {
#'           DATATYPE  H5T_STD_I32LE
#'           DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
#'           DATA {
#'             (0): <position of 1st nt of feature stop codon from orf_gff_file>,
#'                  <position of 2nd nt of feature stop codon from orf_gff_file>,
#'                  <position of 3rd nt of feature stop codon from orf_gff_file>
#'           }
#'         }
#'         DATASET "data" {
#'           DATATYPE  H5T_STD_I32LE
#'           DATASPACE  SIMPLE { ( <sequence_length>, <read_length> ) /
#'                               ( <sequence_length>, <read_length> ) }
#'           DATA {
#'             (0, 0): <see below>
#'             (0, <m>): <see below>
#'             (0, <n>): <see below>
#'             ...
#'             (<p>, 0): <see below>
#'             (<p>, <m>): <see below>
#'             (<p>, <n>): <see below>
#'             ...
#'             (<sequence_length - 1>, 0): <see below>
#'             (<sequence_length - 1>, <m>): <see below>
#'             (<sequence_length - 1>, <n>): <see below>
#'           }
#'         }
#'       }
#'     }
#'   }
#'    ...
#' }
#' ```
#'
#' where:
#'
#' * `gene`: Gene ID from `orf_gff_file` `primary_id` attribute.
#' * `read_length`: `max_read_length` - `min_read_length` + 1
#' * 0 < `m` < `n` < `read_length`. Exact values of `m` and `n` may
#'   differ across specific `ATTRIBUTE` and `DATA` items.
#' * `reads_by_len`:
#'   - `reads_by_len[i]` = number of alignments in `bam_file` which
#'     have `Flag` equal to 0 or 256 and length equal to
#'     lengths[i]`. This equals the sum of `DATA[*, i]` i.e. sum
#'     across all positions for a specific read length.
#' * `sequence_length`: position of final nucleotide of UTR3 from
#'    `orf_gff_file` (equal to length of sequence from `bam_file`
#'    header `LN` value).
#'   - If `is_riboviz_gff` is false when this is equal to `buffer` +
#'     feature length + `buffer`.
#' * `DATASET "data"`:
#'   - 0 <= `p` <= `sequence_length - 1`
#'   - `DATA[p, i]` equals 1 if there is an alignment in `bam_file`
#'     at position `p`+1 which has length equal to `lengths[i]` and
#'     lignment has `Flag` value 0 or 256; 0 otherwise.
#'
#' @export
BamToH5 <- function(bam_file, orf_gff_file, feature, min_read_length,
  max_read_length, buffer, primary_id, secondary_id, dataset,
  stop_in_feature, is_riboviz_gff, hd_file, num_processes) {
  read_lengths <- min_read_length:max_read_length

  # Read in the positions of all exons/genes from GFF and subset
  # feature locations.
  gff <- rtracklayer::readGFFAsGRanges(orf_gff_file)
  gff_features <- gff[gff$type == feature]
  if (length(gff_features) == 0) {
    stop(paste("No", feature, "features found in", orf_gff_file))
  }
  if (!is_riboviz_gff) {
    # GFF does not contain UTR5 or UTR3 elements.
    gff <- gff_features
  }

  # Check primary_id and secondary_id exist.
  gff_col_names <- colnames(rtracklayer::mcols(gff))
  if (!(primary_id %in% gff_col_names)) {
    stop(paste0("primary_id ", primary_id, " is not a GFF attribute"))
  }
  if (!is.na(secondary_id)) {
    if (!(secondary_id %in% gff_col_names)) {
      stop(paste0("Attribute secondary_id ", secondary_id,
           " is not a GFF attribute"))
    }
  }

  # Read in the list of genes from GFF.
  genes <- unique(rtracklayer::mcols(gff)[primary_id][, 1])

  if (!is.na(secondary_id)) {
    alt_genes <- as.list(unique(rtracklayer::mcols(gff)[secondary_id][, 1]))
    names(alt_genes) <- genes
  }
  gff_pid <- rtracklayer::mcols(gff)[primary_id][, 1]

  # Map the reads to individual nucleotide position for each gene.
  print(paste("Mapping reads to feature", feature))
  read_counts <- mclapply(
    genes,
    function(gene) {
      # Select gene location.
      gene_location <- gff[gff_pid == gene]
      # Restrict seqlevels to avoid scanBam error.
      gene_seq_name <- as.character(GenomeInfoDb::seqnames(gene_location)[1])
      GenomeInfoDb::seqlevels(gene_location) <- gene_seq_name
      if (is_riboviz_gff) {
        # GFF contains UTR5, feature, and UTR3 elements.
        buffer_left <- rtracklayer::width(
          gene_location[gene_location$type == "UTR5"])
        buffer_right <- rtracklayer::width(
          gene_location[gene_location$type == "UTR3"])
        gene_location <- gene_location[gene_location$type == feature]
      } else {
        # GFF does not contain UTR5 or UTR3 elements.
        buffer_left <- buffer
        buffer_right <- buffer
      }
      # Get reads to list.
      ReadsToCountMatrix(gene_location = gene_location,
                         bam_file = bam_file,
                         read_lengths = read_lengths,
                         left_flank = buffer_left,
                         right_flank = buffer_right,
                         mult_exon = TRUE)
    },
    mc.cores = num_processes
  )
  names(read_counts) <- genes

  # Save HDF5.
  print(paste("Saving mapped reads in H5 file", hd_file))

  # Create the output HDF5 file.
  rhdf5::h5createFile(hd_file)

  processes <- 1:num_processes
  gene_ids <- seq_len(length(genes))

  # Write base group.
  # Close the file after use to prevent open the same file concurrently.
  fid <- rhdf5::H5Fopen(hd_file)
  base_gid <- rhdf5::H5Gopen(fid, "/")
  rhdf5::H5Fclose(fid)

  for (gene_id in gene_ids) {
    # Create a symbolic link to link HDF5 files back to target HDF5
    # file.
    tmp_hd_file <- paste(hd_file, (gene_id %% num_processes) + 1, sep = ".")
    gene <- genes[gene_id]
    rhdf5::H5Lcreate_external(tmp_hd_file, gene, base_gid, gene)
  }
  for (pid in processes) {
    # Create HDF5 files for each process, so different processes can
    # write into different files. The reason for creating symbolic
    # links first then creating the HDF5 files is that HDF5 will use
    # the absolute path if the source file exists when creating the
    # link. However, we want to use the relative path, so we have to
    # first create link then create file.
    tmp_hd_file <- paste(hd_file, pid, sep = ".")
    rhdf5::h5createFile(tmp_hd_file)
  }
  # Parallelize the processes by writing in different HDF5 files.
  execution_return_list <- mclapply(
    processes,
    function(pid) {
      for (gene_id in gene_ids) {
        if ((gene_id %% num_processes) != pid - 1) {
            # We only store some gene in each process, so skip the
            # irrelevant genes.
          next
        }
        gene <- genes[gene_id]
        # Get matrix of read counts by position and length for the gene.
        gene_read_counts <- read_counts[[gene]]
        # calculate the HDF5 file to write in
        tmp_hd_file <- paste(hd_file, (gene_id %% num_processes) + 1, sep = ".")
        # Get gene location from GFF.
        gene_location <- gff[gff_pid == gene]
        # Get location of start and stop codon nucleotides.
        start_codon_loc <- rtracklayer::start(gene_location[gene_location$type == feature])
        start_cod <- start_codon_loc:(start_codon_loc + 2)
        stop_codon_offset <- 2
        if ((!is_riboviz_gff) && (! stop_in_feature)) {
          # GFF does not contain UTR5 or UTR3 elements.
          stop_codon_offset <- -1
        }
        stop_codon_loc <- rtracklayer::end(gene_location[gene_location$type == feature]) - stop_codon_offset
        stop_cod <- stop_codon_loc:(stop_codon_loc + 2)

        fid <- rhdf5::H5Fopen(tmp_hd_file)
        # Create H5 groups for gene.
        rhdf5::h5createGroup(fid, gene)
        rhdf5::h5createGroup(fid, paste(gene, dataset, sep = "/"))
        rhdf5::h5createGroup(fid, paste(gene, dataset, "reads", sep = "/"))
        mapped_reads <- paste(gene, dataset, "reads", sep = "/")

        gid <- rhdf5::H5Gopen(fid, mapped_reads)
        # Specify then write gene attributes.
        rhdf5::h5createAttribute(gid, "reads_total", c(1, 1))
        rhdf5::h5createAttribute(gid, "buffer_left", c(1, 1))
        rhdf5::h5createAttribute(gid, "buffer_right", c(1, 1))
        rhdf5::h5createAttribute(gid, "start_codon_pos", c(1, 3))
        rhdf5::h5createAttribute(gid, "stop_codon_pos", c(1, 3))
        rhdf5::h5createAttribute(gid, "reads_by_len", c(1, length(read_lengths)))
        rhdf5::h5createAttribute(gid, "lengths", c(1, length(read_lengths)))
        rhdf5:::h5writeAttribute.integer(sum(gene_read_counts),
                                         gid,
                                         name = "reads_total")
        # Though start_codon_loc is an integer, start_codon_loc - 1
        # is a double, so cast back to integer so H5 type is
        # H5T_STD_I32LE and not H5T_IEEE_F64LE.
        rhdf5:::h5writeAttribute.integer(as.integer(start_codon_loc - 1),
                                         gid,
                                         name = "buffer_left")
        rhdf5:::h5writeAttribute.integer((ncol(gene_read_counts) - stop_cod[3]),
                                         gid,
                                         name = "buffer_right")
        rhdf5:::h5writeAttribute.integer(start_cod, gid, name = "start_codon_pos")
        rhdf5:::h5writeAttribute.integer(stop_cod, gid, name = "stop_codon_pos")
        rhdf5:::h5writeAttribute.integer(read_lengths, gid, name = "lengths")
        rhdf5:::h5writeAttribute.integer(apply(gene_read_counts, 1, sum),
                                 gid,
                                 name = "reads_by_len")
        # Specify a dataset within the gene group to store the values and
        # degree of compression, then write the dataset.
        read_data <- paste(mapped_reads, "data", sep = "/")
        rhdf5::h5createDataset(fid,
                        read_data,
                        dim(gene_read_counts),
                        storage.mode = "integer",
                        chunk = c(1, ncol(gene_read_counts)),
                        level = 7)
        rhdf5::h5write(gene_read_counts, fid, name = read_data, start = c(1, 1))
        rhdf5::H5Gclose(gid)
        rhdf5::H5Fclose(fid)
      }
    },
    mc.cores = num_processes,
    # Use preschedule to prevent more than 1 process write the same
    # file at the same time.
    mc.preschedule = TRUE
  )
  # Create symbolic link with alternate IDs, if required.
  if (!is.na(secondary_id)) {
    fid <- rhdf5::H5Fopen(hd_file)
    base_gid <- rhdf5::H5Gopen(fid, "/")
    for (gene in genes) {
      if (alt_genes[[gene]] != gene) {
        rhdf5::H5Lcreate_external(hd_file, gene, base_gid, alt_genes[[gene]])
      }
    }
    rhdf5::H5Gclose(base_gid)
    rhdf5::H5Fclose(fid)
  }
  rhdf5::H5close()
}
