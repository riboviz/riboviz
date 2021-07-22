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

# Define UI ---------------------------------------------------------------

ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  shinyjs::inlineCSS(appCSS),
  
  titlePanel("Riboviz: Generate YAML configuration file"),
  
  p(tags$a(href="https://github.com/riboviz/riboviz",
           "Riboviz GitHub repository")),
  p("For more details about RiboViz workflow configuration parameters: ",
    tags$a(href="https://github.com/riboviz/riboviz/blob/main/docs/user/prep-riboviz-config.md",
           "click here")),
  
  p(),
  p(label_mandatory(""), "Mandatory fields"),
  
  navlistPanel(
    
    "Filenames and locations",
    tabPanel("Sample FASTQ files",
             selectInput("fq_type", "Input FASTQ file(s) is/are:", 
                         c("singleplex", "multiplex")),
             conditionalPanel(condition="input.fq_type == 'singleplex'",
                              helpText("Enter one sample per line as"),
                              helpText("<sample name>: <FASTQ filename>"),
                              helpText("ex. WTnone: SRR1042855_s1mi.fastq.gz"),
                              textAreaInput("fq_files", 
                                            "Sample names and FASTQ filenames (relative to input directory)", "", 
                                            rows=5, width="100%"),
             ),
             conditionalPanel(condition="input.fq_type == 'multiplex'",
                              helpText("Filenames should be relative to input directory"),
                              textInput("multiplex_fq_files", "Multiplexed FASTQ filename", width="100%"),
                              textInput("sample_sheet", "Tab-separated sample sheet filename", width="100%"),
                              helpText("Sample sheet should contain columns 'SampleID' and 'TagRead' (barcode)")
             )
    ),
    tabPanel("Other inputs",
             helpText("Unless otherwise specified, paths should be relative to riboviz directory."),
             textInput("dir_in", label_mandatory("Directory for input files"), width="100%"),
             textInput("dir_index", "Directory for built indices", "index", width="100%"),
             textInput("dir_tmp", "Directory for intermediate files", "tmp", width="100%"),
             textInput("dir_out", "Directory for for output files", "output", width="100%"),
             textInput("asite_disp_length_file", "TSV file with displacement lengths for A-site assignment", width="100%"),
             textInput("features_file", "Features to correlate with ORFs", width="100%"),
             textInput("orf_fasta_file", label_mandatory("FASTA file containing transcript sequences (CDS and flanking regions"), width="100%"),
             textInput("orf_gff_file", label_mandatory("GTF/GFF3 file corresponding to ORF FASTA file"), width="100%"),
             textInput("rrna_fasta_file", label_mandatory("FASTA file containing contaminant (rRNA, etc.) sequences"), width="100%"),
             textInput("codon_positions_file", "Codon positions within each gene (RData file)", width="100%"),
             textInput("t_rna_file", "tRNA estimates (tab-separated values file", width="100%")
    ),
    
    "Riboviz run specifications",
    tabPanel("Riboviz parameters",
             
             helpText("When running Riboviz: "),
             checkboxInput("validate_only", "Only validate configuration without running workflow", F),
             helpText("When validating configuration: "),
             checkboxInput("skip_inputs", "Skip checks for existence of data files", F),
             
             h4("Parameters"),
             textInput("dataset", "Dataset name", "dataset", width="100%"),
             textInput("orf_index_prefix", label_mandatory("Prefix for ORF index files (relative to index directory)"), width="100%"),
             textInput("rrna_index_prefix", label_mandatory("Prefix for contaminant index files (relative to index directory)"), width="100%"),
             textInput("adapters", label_mandatory("Illumina sequencing adapter to remove"), width="100%"),
             textInput("umi_regexp", "Regular expression to extract barcodes and UMIs using UMI-tools", width="100%"),
             numericInput("count_threshold", "Minimum read count to keep gene in report", 1, width="100%"),
             numericInput("min_read_length", "Minimum read length in H5 output", 10, width="100%"),
             numericInput("max_read_length", "Maximum read length in H5 output", 50, width="100%"),
             textInput("primary_id", "Primary gene IDs (ex. YAL001C, YAL003W)", "Name", width="100%"), # ???
             numericInput("num_processes", "Number of processes to parallelize over", 1, width="100%"),
             textInput("samsort_memory", "Memory to allocate to samtools", "null", width="100%"),
             numericInput("buffer", "Length of flanking region around CDS", 250), width="100%",
             textInput("secondary_id", "Secondary gene IDs (ex. COX1, EFB1)", "NULL", width="100%") # ???
    ),
    tabPanel("Other options",
             checkboxInput("rpf", "Dataset is an RPF (vs. mRNA) datset", F, width="100%"), # ???
             textInput("feature", "Feature type", "CDS", width="100%"), # ???
             # checkboxInput("stop_in_cds", "Stop codons are part of CDS annotations in GFF/GTF3", F),
             checkboxInput("stop_in_feature", "Stop codons are part of feature annotations in GFF/GTF3", F, width="100%"),
             checkboxInput("build_indices", "Build indices for aligner", T, width="100%"),
             checkboxInput("count_reads", "Count reads in input, temporary, and output files", T, width="100%"),
             checkboxInput("trim_5p_mismatches", "Trim mismatched 5' base", T, width="100%"),
             checkboxInput("dedup_umis", "Deduplicate reads using UMI-tools", F, width="100%"),
             checkboxInput("extract_umis", "Extract UMIs after adapter trimming", F, width="100%"),
             checkboxInput("group_umis", "Summarize UMI groups pre- and post-deduplication", F, width="100%"),
             checkboxInput("dedup_stats", "Output UMI deduplication statistics", T, width="100%"),
             checkboxInput("do_pos_sp_nt_freq", "Calculate position-specific nucleotide frequency", T, width="100%"),
             checkboxInput("is_riboviz_gff", "GFF file contains 3 elements (UTR5, CDS, UTR3) per gene", T, width="100%"),
             checkboxInput("make_bedgraph", "Output bedgraph data files in addition to H5 files", T, width="100%"),
             checkboxInput("output_pdfs", "Generate .pdfs for sample-related plots", T, width="100%"),
             checkboxInput("run_static_html", "Run static HTML visualization per sample", T, width="100%"),
             checkboxInput("publish_index_tmp", "Copy index and temporary files from Nextflow's <work/> directory", F, width="100%"),
    ),
    
    "Nextflow and batch job specifications",
    tabPanel("Nextflow options",
             p("For more details: ", 
               tags$a(href="https://github.com/riboviz/riboviz/blob/main/docs/user/create-job-script.md",
                      "click here"),),
             textInput("nextflow_work_dir", "Nextflow work directory", "work", width="100%"),
             textInput("nextflow_dag_file", "Nextflow DAG filename", "nextflow-dag.html", width="100%"),
             textInput("nextflow_report_file", "Nextflow report file", "nextflow-report.html", width="100%"),
             textInput("nextflow_timeline_file", "Nextflow timeline filename", "nextflow-timeline.html", width="100%"),
             textInput("nextflow_trace_file", "Nextflow trace filename", "nextflow-trace.tsv", width="100%")
    ),
    tabPanel("Batch job parameters",
             p("For more details: ", 
               tags$a(href="https://github.com/riboviz/riboviz/blob/main/docs/user/create-job-script.md",
                      "click here"),),
             textInput("job_name", "Name of batch job", "riboviz", width="100%"),
             textInput("job_email", "Email address for batch job events", width="100%"),
             textInput("job_email_events", "Events to trigger emails", "beas", width="100%"),
             textInput("job_memory", "Requested memory", "8G", width="100%"),
             numericInput("job_num_cpus", "Requested number of CPUs", "4", width="100%"),
             textInput("job_parallel_env", "Requested parallel environment", "mpi", width="100%"),
             textInput("job_runtime", "Maximum runtime", "48:00:00", width="100%")
    ),
    
    "-----",
    
    tabPanel("Download YAML",
             helpText("Copy text below or click download button to generate new yaml config file."),
             downloadButton("download_yaml", label="Download"),
             p(),
             verbatimTextOutput("text")
    ),
    
    widths=c(4,8)
  )
)

# Define server logic -----------------------------------------------------

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
    output_text <- reactive({
      yaml_params <- sort(names(input))
      yaml_params <- subset(yaml_params, yaml_params != "fq_type")
      sapply(yaml_params,
             function(param) {
               if(param=="fq_files") {
                 if(input$fq_files == "") { # no samples provided
                   "fq_files:"
                 } else {
                   paste0("fq_files: \n  ",
                          gsub("\n", "\n  ",
                               input$fq_files))
                 }
               } else {
                 paste0(param, ": ", input[[param]])
               }
             })
    })
    
    # output yaml text
    output$text <- renderPrint({ cat(output_text(), sep="\n") })
    
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
}

# Run app -----------------------------------------------------------------

shinyApp(ui=ui, server=server)