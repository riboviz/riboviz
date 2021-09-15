source(here::here("rscripts", "read_count_functions.R"))
source(here::here("rscripts", "stats_figs_block_functions.R"))

suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))

#### Create base table functions ####

#' FilterForFrameFunction
#' FilterForFrameFunction() is an alternate way of mapping reads to codons.
#' Unlike SnapToCodon, FilterForFrameFunction only returns reads mapping to
#' the first nucleotide of a codon for the desired frame.
#'
#' @param gene - Gene name to pull out read counts for
#' @param gff_df - The data frame version of the GFF3 file
#' @param dataset - The name of dataset stored in .h5 file.
#' @param hd_file - The path to the .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param min_read_length - An integer, the minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param asite_displacement_length - The length from the end of a read to the A site
#' @param filtering_frame -  The frame from which to select reads from
#'
#' @return - a list of numeric values (read counts) for a reading frame (0/1/2)
#'
#' @examples FilterForFrameFunction(gene, dataset, hd_file,
#'                                  min_read_length,
#'                                  asite_displacement_length,
#'                                  gff_df,
#'                                  filtering_frame)

FilterForFrameFunction <- function(gene, gff_df, dataset, hd_file,
                                   min_read_length, asite_displacement_length,
                                   filtering_frame) {
  
  # Get the matrix of read counts
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hd_file)
  
  # Assign reads to their A site nt, determined by read length
  reads_asitepos <- CalcAsiteFixed(reads_pos_length,  
                                   min_read_length,
                                   asite_displacement_length
  )
  
  # Get the gff rows for the gene being studied
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene)
  
  # Get the position of the start codon
  left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                   type == "CDS") %>%  select(start))
  
  # Get the position of the stop codon
  right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                    type == "CDS") %>%  select(end))
  
  # Extract the reads that map to the CDS
  cds_reads <- reads_asitepos[left:right]
  # ie: num [1:621] 811 460 2978 429 251 ...
  
  cds_length <- length(cds_reads)/3
  
  # Create a tibble, assigning a frame to each nt, so the first nt in each frame
  # has the corresponding frame identity
  cds_reads_frames <- tibble(Count = cds_reads,
                             Frame = rep(c(0, 1, 2), times = cds_length)
  )
  # Example
  # > str(cds_frames)
  # tibble [621 x 2] (S3: tbl_df/tbl/data.frame)
  # $ Count: num [1:621] 811 460 2978 429 251 ...
  # $ Frame: num [1:621] 0 1 2 0 1 2 0 1 2 0 ...
  
  filtered_for_frame <- dplyr::filter(cds_reads_frames, 
                                      Frame == filtering_frame)
  # > str(cds_frames)
  # tibble [621 x 2] (S3: tbl_df/tbl/data.frame)
  # $ Count: num [1:621] 811 460 2978 429 251 ...
  # $ Frame: num [1:621] 0 1 2 0 1 2 0 1 2 0 ...
  
  filtered_counts <- filtered_for_frame$Count
  # > str(filtered_counts)
  # num [1:207] 811 429 488 102 994 146 173 762 13 176 ...
  
  return(filtered_counts)
  
}



#' GetAllCodonPosCounts1Gene()
#' GetAllCodonPosCounts1Gene works for one gene, and calculates the read
#' counts at each codon position in that gene.
#' @param gene 
#' @param gff_df 
#' @param asite_disp_path 
#' @param dataset 
#' @param hd_file 
#' @param min_read_length 
#' @param filter_for_frame 
#' @param filtering_frame 
#' @param snapdisp 
#'
#' @return
#' @export
#'
#' @examples
GetAllCodonPosCounts1Gene <- function(gene, gff_df, asite_disp_path, dataset, hd_file, min_read_length,
                                      filter_for_frame, filtering_frame, snapdisp) {
  
  # Get the gff rows for the gene being studied
  subset_gff_df_by_gene <- dplyr::filter(.data = gff_df, seqnames == gene)
  
  # Get the position of the start codon
  left <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                   type == "CDS") %>%  select(start))
  
  # Get the position of the stop codon
  right <- as.numeric(dplyr::filter(.data = subset_gff_df_by_gene,
                                    type == "CDS") %>%  select(end))
  
  # Get the Asite displacement for reads of different lengths
  asite_displacement_length <- suppressMessages(ReadAsiteDisplacementLengthFromFile(asite_disp_path))
  
  # Calculate reads at codon positions in desired way, decided by
  # filter_for_frame being TRUE or FAlSE (default = TRUE)
  
  if(filter_for_frame == FALSE){
    
    # Use SnapToCodon to assign reads to positions
    # This adds the reads from all three nucleotides in the codon
    codon_counts_1_gene <- GetGeneCodonPosReads1dsnap(gene, dataset, hd_file,
                                                      min_read_length,
                                                      asite_displacement_length,
                                                      left, right,
                                                      snapdisp)
    # > str(codon_counts_1_gene)
    # num [1:207] 4249 825 1017 1176 1116 ...
    
  } else {
    
    # Use FilterForFrameFunction to assign reads to positions
    # This uses only the reads from the nucleotide in the codon of a designated frame
    codon_counts_1_gene <- FilterForFrameFunction(gene, gff_df, dataset, hd_file,
                                                  min_read_length,
                                                  asite_displacement_length,
                                                  filtering_frame)
    # str(codon_counts_1_gene)
    # num [1:207] 811 429 488 102 994 146 173 762 13 176
    
  }
  
  # Add codon positions
  codon_pos_counts <- tibble(Gene = gene,
                             PosCodon = 1:length(codon_counts_1_gene),
                             Count = codon_counts_1_gene)
  
  as.data.frame(codon_pos_counts, row.names = NULL, optional = FALSE)
  
  return(codon_pos_counts)
  
}




#' GetAllCodonPosCounts
#' 
#' Applies the GetGeneCodonPosReads1dsnap() function to a list of genes and
#' generates a tidy data frame (tibble) which contains the counts for all genes
#'
#' @param gene_names - The list of all the genes in the sample
#' @param gff_df - The data frame version of the GFF3 file
#' @param dataset - The name of dataset stored in .h5 file.
#' @param hd_file - The path to the .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param min_read_length - An integer, the minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param filter_for_frame - Decide which method to use for assigning reads, TRUE or FALSE
#' @param filtering_frame - The frame from which to select reads from
#' @param snapdisp - frame to filter to when using SnapToCodon
#'
#' @return - A tibble listing the Gene, PosCodon and Count
#' @export
#'
#' @examples total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, gff_df, 
#' dataset, hd_file, min_read_length, filter_for_frame, filtering_frame,  snapdisp)
#' 
GetAllCodonPosCounts <- function(gene_names, gff_df, asite_disp_path, dataset, hd_file, min_read_length,
                                 filter_for_frame, filtering_frame, snapdisp) {
  
  # purrr::map is used to apply GetAllCodonPosCounts1Gene to all genes in the
  # sample
  total_codon_pos_counts <- purrr::map_dfr(.x = gene_names,
                                           .f = GetAllCodonPosCounts1Gene,
                                           gff_df,
                                           asite_disp_path,
                                           dataset,
                                           hd_file,
                                           min_read_length,
                                           filtering_frame = filtering_frame,
                                           filter_for_frame = filter_for_frame,
                                           snapdisp = snapdisp
  )
  
  return(total_codon_pos_counts)
}



#' AddCodonNamesToCodonPosCounts 
#' AddCodonNamesToCodonPosCounts() Adds codon identities to position and counts
#'
#' @param gene_names - The list of all the genes in the sample
#' @param gff_df - The gff file in a dataframe format
#' @param dataset - The name of dataset stored in .h5 file.
#' @param hd_file - The path to the .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param min_read_length -  An integer, the minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param filter_for_frame - Decides which method of asite assignment to use
#' @param filtering_frame - The frame from which to select reads from
#' @param snapdisp - An integer for additional displacement in the snapping
#' @param yeast_codon_pos_i200 - A list of all codons and positions from the sample being studied
#'
#' @return - a tibble which contains the columns "Gene", "PosCodon", "Count",
#' and "codon" for a list of genes
#' @export
#'
#' @examples AddCodonNamesToCodonPosCounts(gene_names, gff_df, dataset, hd_file,
#' min_read_length, filter_for_frame,filtering_frame, snapdisp,yeast_codon_pos_i200)
#' 
AddCodonNamesToCodonPosCounts <- function(gene_names, gff_df, asite_disp_path, dataset, hd_file, 
                                          min_read_length, filter_for_frame,
                                          filtering_frame, snapdisp, yeast_codon_pos_i200) {
  
  # Create tibble of read counts and positions
  total_codon_pos_counts <- GetAllCodonPosCounts(gene_names, gff_df, asite_disp_path, dataset, 
                                                 hd_file, min_read_length,
                                                 filter_for_frame, 
                                                 filtering_frame, snapdisp)
  
  # Add codon identities
  transcript_tibbles <- left_join(total_codon_pos_counts, yeast_codon_pos_i200,
                                  by = c("PosCodon", "Gene"), keep = FALSE,
                                  copy = TRUE)
  
  # turn into a tibble
  transcript_gene_poscodon_frame <- tibble(
    Gene = transcript_tibbles$Gene,
    PosCodon = transcript_tibbles$PosCodon,
    Count = transcript_tibbles$Count,
    Codon = transcript_tibbles$Codon
  )
  
  return(transcript_gene_poscodon_frame)
}



#### Slice window around interesting features #####



#' TranscriptForOneGene
#' Take an individual gene as an input, then filter for the codon of interest
#'  on the gene being investigated
#' @param gene 
#' @param dataset 
#' @param hd_file 
#' @param min_read_length 
#' @param gff_df 
#' @param yeast_codon_pos_i200 
#' @param feature_of_interest 
#'
#' @return
#' @export
#'
#' @examples
TranscriptForOneGene <- function(gene, dataset, hd_file, 
                                 min_read_length,   
                                 gff_df, yeast_codon_pos_i200, 
                                 feature_of_interest) {
  
  # Filter transcript_gene_poscodon_frame to only include features of interest
  interesting_feature_tibble <- dplyr::filter(
    transcript_gene_poscodon_frame, 
    Codon == feature_of_interest)
  
  # select only the gene that the function is being applied to in this 
  # purrr::map cycle
  transcript_for_one_gene <- dplyr::filter(
    interesting_feature_tibble, Gene == gene)
  
  return(transcript_for_one_gene)
  
}


#' ExpandRegions
#' For each occurrence of the feature_of_interest on the gene being studied, 
#' slice out a window of positions around the feature_of_interest from 
#'  transcript_gene_poscodon_frame 
#'  return an empty tibble if the desired region hangs over the edge of the 
#'  coding region
#'  This if statement ensures that feature tibbles containing more or less 
#'  positions than expected are discarded
#' @param transcript_for_one_gene 
#' @param transcript_gene_poscodon_frame 
#' @param gene 
#' @param dataset 
#' @param hd_file 
#' @param expand_width 
#'
#' @return
#' @export
#'
#' @examples
ExpandRegions <- function(transcript_for_one_gene, 
                          transcript_gene_poscodon_frame, gene, dataset,
                          hd_file, expand_width){
  
  # save applied transcript to another object remove problems with purrr::map
  interesting_features <- transcript_for_one_gene
  
  # Filter transcript_gene_pos_poscodon_gene_interest to just the gene being
  # used in this purrr::map cycle
  transcript_gene_pos_poscodon_gene_interest <- dplyr::filter(
    transcript_gene_poscodon_frame, Gene == gene)
  
  # get the length of the gene being studied
  gene_length <- filter(gff_df,
                        gff_df$type == "CDS" & gff_df$Name == gene)$width
  
  # if the position of the feature_of_interest is near the start or stop codon
  # return a NULL in resulting list of tibbles.
  if (interesting_features <= expand_width | interesting_features +
      expand_width > gene_length/3) {
    
    return()
    
  } else {
    
    # Slice out a window of positions around the feature_of_interest
    output_feature_info <- tibble(
      dplyr::slice(transcript_gene_pos_poscodon_gene_interest,
                   (interesting_features - expand_width):(interesting_features + expand_width), 
                   each = FALSE),
      Rel_Pos =  seq(- expand_width, expand_width)
    )
    
    if(dim(output_feature_info)[1] == (2*expand_width + 1)){
      
      return(output_feature_info)
    }else{
      return()
    }
  }
}


#' AllGeneInterestingFeatures
#'
#' @param gene 
#' @param gene_names 
#' @param dataset 
#' @param hd_file 
#' @param min_read_length 
#' @param gff_df 
#' @param yeast_codon_pos_i200 
#' @param feature_of_interest 
#' @param transcript_gene_poscodon_frame 
#'
#' @return
#' @export
#'
#' @examples
AllGeneInterestingFeatures <- function(gene, gene_names, dataset, hd_file, 
                                       min_read_length, gff_df,
                                       yeast_codon_pos_i200, feature_of_interest, 
                                       transcript_gene_poscodon_frame) {
  
  
  # select features of interest occurring on the gene being that the function 
  # is being applied to in this purrr::map cycle 
  transcript_for_one_gene <- TranscriptForOneGene(gene, dataset, hd_file, 
                                                  min_read_length, gff_df, 
                                                  yeast_codon_pos_i200,
                                                  feature_of_interest)
  
  
  # Use transcript_for_one_gene as an input to ExpandRegions
  
  # Apply ExpandRegion() to each occurrence of the feature_of_interest 
  # on the gene being studied.
  output_feature_info <- purrr::map(.x = transcript_for_one_gene$PosCodon, 
                                    .f = ExpandRegions, 
                                    transcript_gene_poscodon_frame,
                                    gene,
                                    dataset,
                                    hd_file,
                                    expand_width)
  
  
  return(output_feature_info)
}




#' ExpandFeatureRegionAllGenes
#' this function takes a slice out of
#' transcript_gene_poscodon_frame, centered on the position of
#' the feature_of_interest
#' It then sets the position of the feature_of_interest to 0, and changes the
#' positions of adjacent codons to be relative to the feature_of_interest
#' This is done for all occurrences of the feature_of_interest on all genes in
#' the sample provided.
#'
#' @param gene_names - The list of all the genes in the sample
#' @param gff_df - The gff file in a dataframe format
#' @param dataset - The name of dataset stored in .h5 file.
#' @param hd_file - The path to the .h5 hdf5 file holding read data for all genes,
#' created from BAM files for dataset samples.
#' @param min_read_length - An integer, the minimum read length in H5 output;
#' Default = 10 (set in generate_stats_figs.R from yaml)
#' @param filter_for_frame - Decides which method of asite assignment to use
#' @param filtering_frame - The frame from which to select reads from
#' @param snapdisp - An integer for additional displacement in the snapping
#' @param yeast_codon_pos_i200 - A list of all codons and positions from the sample being studied
#' @param feature_of_interest - The feature being studied
#' @param expand_width - the number of codons to take a slice of eithr side of 
#' occurrences of the feature_of_interest
#'
#' @return - A list of tibbles, each containing the columns "Gene", "PosCodon", 
#' "Count", "Codon, "Rel_Pos" for a window of positions determined by the 
#' expand_width around an occurance of the feature_of_interest 
#' @export
#'
#' @examples     ExpandFeatureRegionAllGenes(gene_names, gff_df, dataset, hd_file, 
#' min_read_length, filter_for_frame, filtering_frame, snapdisp, 
#' yeast_codon_pos_i200,feature_of_interest,expand_width)
ExpandFeatureRegionAllGenes <- function(gene_names, gff_df, asite_disp_path,
                                        dataset, hd_file, 
                                        min_read_length,   
                                        filter_for_frame,
                                        filtering_frame, snapdisp, yeast_codon_pos_i200,
                                        feature_of_interest,
                                        expand_width) {
  
  # generate transcript_gene_poscodon_frame, tests and descriptions above
  transcript_gene_poscodon_frame <- AddCodonNamesToCodonPosCounts(gene_names,
                                                                  gff_df,
                                                                  asite_disp_path,
                                                                  dataset, 
                                                                  hd_file, 
                                                                  min_read_length, 
                                                                  filter_for_frame,
                                                                  filtering_frame, 
                                                                  snapdisp,
                                                                  yeast_codon_pos_i200)
  
  
  # Apply AllGeneInterestingFeatures to each gene in the sample
  output_feature_info <- purrr::map(.x = gene_names, 
                                    .f = AllGeneInterestingFeatures,
                                    gene_names,
                                    dataset,
                                    hd_file, 
                                    min_read_length,
                                    gff_df,
                                    yeast_codon_pos_i200,
                                    feature_of_interest,
                                    transcript_gene_poscodon_frame)
  
  # produces a list for each gene, containing a list for each occurrence of the 
  # feature_of_interest
  
  # Unlist to produce one list, containing each occurrence of the feature_of_interest
  output_feature_info <- unlist(output_feature_info, recursive = F)
  
  # remove NULLS, which represent the feature_of_interest occurring within one 
  # expand_width of the UTRs
  output_feature_info <- output_feature_info[!sapply(output_feature_info, is.null)]
  
  return(output_feature_info)
}

#### Functions for normalization #####


#' SetNaNToZero
#'
#' @param relcount_values 
#'
#' @return
#' @export
#'
#' @examples
SetNaNToZero <- function(relcount_values){
  
  if(is.nan(relcount_values)){
    relcount_values <- 0
  }else{
    relcount_values <- relcount_values 
  }
}


#' CheckForNaN
#'
#' @param normalized_expand_tibble 
#'
#' @return
#' @export
#'
#' @examples
CheckForNaN <- function(normalized_expand_tibble){
  
  relcount_values <- unlist(normalized_expand_tibble$RelCount)
  
  relcount_values <- unlist(purrr::map(.x = relcount_values, .f = SetNaNToZero))
  
  dplyr::mutate(.x, RelCount = relcount_values)
}



#' ExpandedRegionNormalization(): carries out normalization within each expanded frame 
#' Normalization carried out within each expanded frame so that they are comparable 
#' Normalizes the expand_feature_region list generating a RelCount column with the
#' normalization values
#' Normalizes the ExpandFeatureRegion list generating a RelCount column with the 
#' normalization values.
#' As this function is not looped it will only generate one normalized tibble 
#' for each occurrence of the feature_of_interest 
#' 
#' @param .x - The list of tidy format data frames (tibbles) generated 
#' by the function ExpandFeatureRegion
#' @param expand_width - an integer which provides the number of positions on each 
#' side of the feature_of_interest to include in the window
#' 
#' @return the list of tibbles which contain the normalized counts within the 
#' window so that the feature is comparable despite overall varying levels of 
#' expression between genes 
#' @export
#' @example ExpandedRegionNormalization(.x, expand_width)
#' 
ExpandedRegionNormalization <- function(.x, expand_width){
  
  # Normalize the read counts 
  normalized_expand_tibble <- dplyr::mutate(.x, RelCount = Count / sum(Count) * (2 * expand_width + 1))
  
  # Having NaNs in any table will cause problems during the overlay step
  # Changing them to 0 removes this problem, and takes into account that the 
  # NaNs are caused by no reads occurring
  
  CheckForNaN(normalized_expand_tibble)
}


#### Overlay functions ####



#' OverlayedTable(): overlays tidy format data frames (tibbles) to create a single overlayed tibble.
#' 
#' FIXME: Takes normalized_expand_list as its input.
#' 
#' @param normalized_expand_list the output from the looped function ExpandedRegionNormalization()
#' @param expand_width integer which provides the number of positions on each side 
#' of the feature_of_interest to include in the window
#' 
#' @return a tibble which contains the mean counts for each position from the normalized tibbles. 
#' @export
#' 
#' @example  normalized_expand_list <- purrr::map(.x = expand_feature_region,
#'                                      .f = ExpandedRegionNormalization,
#'                                      expand_width)
#' 
#' OverlayedTable(normalized_expand_list, expand_width)

OverlayedTable <- function(normalized_expand_list, expand_width){
  
  # The number of objects inside normalized_expand_list
  number_of_objects <- length(normalized_expand_list)
  
  # Reduces normalized_expand_list to the columns Rel_Pos and RelCount
  result <- lapply(normalized_expand_list, "[", c("Rel_Pos", "RelCount"))
  
  joined_result <- result %>% purrr::reduce(full_join, by = c("Rel_Pos"), sum("RelCount"))
  
  # rejoin values
  joined_rows = joined_result %>% 
    mutate(SumRows = rowSums(select(joined_result, -"Rel_Pos")) / number_of_objects)
  
  # recreate final tibble 
  overlayed_tibbles <- tibble::tibble(
    Rel_Pos = seq(- expand_width, expand_width),
    RelCount = joined_rows$SumRows
  )
}


