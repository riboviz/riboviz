library(shiny)
library(shinyjs)

mandatory_fields <- c("adapters",
                      "dir_in",
                      "orf_fasta_file",
                      "orf_gff_file",
                      "orf_index_prefix",
                      "rrna_fasta_file",
                      "rrna_index_prefix")

label_mandatory <- function(label) {
  tagList(label, span("*", class="mandatory_star"))
}

appCSS <- ".mandatory_star {color: red; }"

server <- function(input, output) {
  observe({
    # check if all mandatory fields have a value
    mandatory_filled <- vapply(mandatory_fields,
                               function(x) {
                                 !is.null(input[[x]]) && 
                                   input[[x]] != ""
                               },
                               logical(1))
    all_mandatory_filled <- all(mandatory_filled)
    
    # enable/disable the download button
    shinyjs::toggleState(id="download_yaml", condition=all_mandatory_filled)
    
    # generate yaml text
    output_text <- c(paste("adapters:", input$adapters,
                           "# Illumina sequencing adapter(s) to remove"),
                     paste("asite_disp_length_file:", input$asite_disp_length_file,
                           "# Table of fixed A-site positions by read length"),
                     paste("buffer:", input$buffer,
                           "# Length of flanking region around the CDS"),
                     paste("build_indices:", as.character(input$build_indices),
                           "# Build indices for aligner? if TRUE, remake indices",
                           "from fasta files"),
                     paste("codon_positions_file:", input$codon_positions_file,
                           "# Codon positions in each gene"),
                     paste("count_reads:", as.character(input$count_reads),
                           "# Scan input, temporary and output files and produce",
                           "counts of reads in each FASTQ, SAM, and BAM file processed?"),
                     paste("count_threshold:", input$count_threshold,
                           "# Remove genes with a read count below this threshold, when",
                           "generating statistics and figures"),
                     paste("dataset:", input$dataset,
                           "# Dataset name"),
                     paste("dedup_stats:", as.character(input$dedup_stats),
                           "# Output UMI deduplication statistics?"),
                     paste("dedup_umis:", as.character(input$dedup_umis),
                           "# Extract UMIs and deduplicate reads if TRUE"),
                     paste("dir_index:", input$dir_index,
                           "# Input directory"),
                     paste("dir_out:", input$dir_out,
                           "# Output directory"),
                     paste("dir_tmp:", input$dir_tmp,
                           "# Intermediate files directory"),
                     paste("do_pos_sp_nt_freq:", as.character(input$do_pos_sp_nt_freq),
                           "# Calculate position-specific nucleotide frequency?"),
                     paste("extract_umis:", as.character(input$extract_umis),
                           "# Extract UMIs if TRUE"),
                     paste("feature:", input$feature,
                           "# Feature type"),
                     paste("features_file:", input$features_file,
                           "# Features to correlate with ORFs"),
                     "fq_files: # fastq files to be processed, relative to dir_in",
                     {
                       if(input$sample_names == "") {
                         ""
                       } else {
                         paste("  ", unlist(strsplit(tmp, split="\n")))
                       }
                     },
                     paste("group_umis:", as.character(input$group_umis),
                           "# Summarise UMI groups before and after deduplication, if TRUE"),
                     paste("is_riboviz_gff:", as.character(input$is_riboviz_gff),
                           "# Does the GFF file contain 3 elements per gene - UTR5, CDS, and UTR3"),
                     paste("job_email_events:", input$job_email_events,
                           "# Events triggering emails about batch job (job submission).",
                           "Any combination of b - begin, e - end, a - abort, s - suspend."),
                     paste("job_email:", input$job_email,
                           "# E-mail address for batch job events (job submission)."),
                     paste("job_memory:", input$job_memory,
                           "# Requested memory for batch job (job submission)."),
                     paste("job_name:", input$job_name,
                           "# Name of batch job (job submission)."),
                     paste("job_num_cpus:", input$job_num_cpus,
                           "# Requested number of CPUs for batch job (job submission)."),
                     paste("job_parallel_env:", input$job_parallel_env,
                           "# Requested parallel environment for batch job (Grid Engine job submission)."),
                     paste("job_runtime:", input$job_runtime,
                           "# Maximum runtime for batch job (job submission)."),
                     paste("make_bedgraph:", input$make_bedgraph,
                           "# Output bedgraph files, as TSV, in addition to h5?"),
                     paste("max_read_length:", input$max_read_length,
                           "# Maximum read length in H5 output"),
                     paste("min_read_length:", input$min_read_length,
                           "# Minimum read lenght in H5 output"),
                     paste("multiplex_fq_files:", input$multiplex_fq_files,
                           "# Multiplexed fastq files to be processed, relative to dir_in"),
                     paste("nextflow_dag_file:", input$nextflow_dag_file,
                           "# Nextflow DAG file (job submission)."),
                     paste("nextflow_report_file:", input$nextflow_report_file,
                           "# Nextflow report file (job submission)."),
                     paste("nextflow_timeline_file:", input$nextflow_timeline_file,
                           "# Nextflow timeline file (job submission)."),
                     paste("nextflow_trace_file:", input$nextflow_trace_file,
                           "# Netflow trace file (job submission)."),
                     paste("nextflow_work_dir:", input$nextflow_work_dir,
                           "# Nextflow work directory (job submission)."),
                     paste("num_processes:", input$num_processes,
                           "# Number of processes to parallelize over"),
                     paste("orf_fasta_file:", input$orf_fasta_file,
                           "# ORF file to align to"),
                     paste("orf_gff_file:", input$orf_gff_file,
                           "# GFF2/GFF3 file for ORFs"),
                     paste("orf_index_prefix:", input$orf_index_prefix,
                           "# ORF idnex file prefix, relative to dir_index"),
                     paste("output_pdfs:", as.character(input$output_pdfs),
                           "# generate .pdfs for sample-related plots"),
                     paste("primary_id:", input$primary_id,
                           "# Primary gene IDs to access the data (YAL001C, YAL003W, etc.)"),
                     paste("publish_index_tmp:", as.character(input$publish_index_tmp),
                           "# Publish index and temporary files to dir_index and dir_tmp? If FALSE, use symlinks."),
                     paste("rpf:", as.character(input$rpf),
                           "# Is the dataset an RPF or mRNA dataset?"),
                     paste("rrna_fasta_file:", input$rrna_fasta_file,
                           "# rRNA file to avoid aligning to"),
                     paste("rrna_index_prefix:", input$rrna_index_prefix,
                           "# rRNA index file prefix, relative to dir_index"),
                     # paste("run_static_html:", as.character(input$run_static_html),
                     #       "# Create static html visualization per sample?"),
                     ### ^ still in vignette config, but not listed on
                     ### https://github.com/riboviz/riboviz/blob/main/docs/user/prep-riboviz-config.md
                     paste("sample_sheet:", input$sample_sheet,
                           "# Sample sheet, TSV file with, at least,",
                           "SampleID and TagRead (barcode) columns"),
                     paste("samsort_memory:", input$samsort_memory,
                           "# Memory to give to 'samtools sort'"),
                     paste("secondary_id:", input$secondary_id,
                           "# Secondary gene IDs to access the data (COX1, EFB1, etc.)"),
                     paste("stop_in_feature:", as.character(input$stop_in_feature),
                           "# Are stop codons part of the feature annotations in GFF?"),
                     paste("trim_5p_mismatches:", as.character(input$trim_5p_mismatches),
                           "# Trim mismatched 5' base?"),
                     paste("t_rna_file:", input$t_rna_file,
                           "# tRNA estimates"),
                     paste("umi_regexp:", input$umi_regexp,
                           "# UMI-tools-compliant regular expression to extract UMIs"),
                     paste("validate_only:", as.character(input$validate_only))
    )
    output_text <- output_text[output_text != ""]
    
    # output yaml text
    output$text <- renderPrint({ cat(output_text, sep="\n") })
    
    # download yaml text
    output$download_yaml <- downloadHandler(
      filename = function() {
        "riboviz_config.yaml"
      },
      content = function(file) {
        writeLines(sub("\n", "", output_text), file)
      }
    )
  })
