#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

/*
===================================
RiboViz ribosome profiling workflow
===================================
*/

def helpMessage() {
    log.info """

    Usage
    -----

        nextflow run prep_riboviz.nf -params-file <CONFIG>.yaml [--help]

    where '<CONFIG>' is a YAML configuration file. The YAML
    configuration parameters are as follows (all are mandatory unless
    stated).

    Configuration parameters can also be provided via the
    command-line in the form '--<PARAMETER>=<VALUE>' (for example
    '--make_bedgraph=FALSE').

    '--help' displays this help information and exits.

    Configuration
    -------------

    Organism data:

    * 'orf_fasta_file': Transcript sequences file containing both
      coding regions and flanking regions (FASTA file)
    * 'orf_gff_file': Matched genome feature file, specifying coding
      sequences locations (start and stop coordinates) within the
      transcripts (GTF/GFF3 file)
    * 'rrna_fasta_file': Ribosomal rRNA and other contaminant
      sequences to avoid aligning to (FASTA file)

    Ribosome profiling data:

    * 'dir_in': Input directory.
    * 'download' Optional dictionary of SRA accession numbers to
      download data directly from SRA. Downloaded data will be placed in 
      directory specified as the input directory. Each item consists of a sample
      name with a filenmae value (e.g. 'WT3AT: SRR1042864'). 
      In the case of multiplexed data, the choice of sample name is arbitrary,
      as the true sample names will be taken from the 'sample_sheet'.
      Currently, user must still provide 'fq_files' or 'multiplex_fq_files'
      as specified, which can be simply named as <ACCESSION>.fastq.gz. Note that
      code has only been tested for accessions starting with the run accession number
      starting with SRR. 
    * Either:
      - 'fq_files': Dictionary of FASTQ files to be processed,
        relative to '<dir_in>'. Each item consists of a sample name
        with a file name value
        (e.g. 'WT3AT: SRR1042864_s1mi.fastq.gz')
    * Or:
      - 'multiplex_fq_files': List with a multiplexed FASTQ file,
        relative to '<dir_in>'. If this is provided then the
        'fq_files' parameter must not be present in the configuration
        and the 'sample_sheet' parameter must be present.
      - 'sample_sheet': A sample sheet, relative to '<dir_in>',
        mandatory if 'multiplex_fq_files' is used (tab-separated
        values file with, at least, 'SampleID' and 'TagRead' (barcode)
        columns)
    * If neither or both of 'fq_files' and 'multiplex_fq_files'
      parameters are provided then the workflow will exit.

    Indexing:

    * 'build_indices': Rebuild indices from FASTA files? (default
      'TRUE'). If 'FALSE' then 'dir_index' is expected to contain the
      index files.
    * 'orf_index_prefix': Prefix for ORF index files, relative to
      '<dir_index>' .
    * 'rrna_index_prefix': Prefix for rRNA index files, relative to
      '<dir_index>'.
    * 'dir_index': Directory to write indexed files to (default 'index')

    Outputs:

    * 'dir_tmp': Directory to write temporary files to (default 'tmp')
    * 'dir_out': Directory to write temporary files to
      (default 'output')

    Adapter trimming:

    * 'adapters': Illumina sequencing adapter(s) to remove.

    Barcode and UMI extraction, deduplication, demultiplexing:

    * 'extract_umis': Extract UMIs after adapter trimming? (default
      'FALSE')
    * 'umi_regexp': UMI-tools-compliant regular expression to extract
      barcodes and UMIs. For details on the regular expression format,
      see UMI-tools documentation on Barcode extraction
      https://umi-tools.readthedocs.io/en/latest/reference/extract.html#barcode-extraction.
      Only required if 'extract_umis' is 'TRUE'.
      - If 'fq_files' are provided then 'umi_regexp' should extract
        only UMIs (i.e. it should contain '<umi>' elements only).
      - If 'multiplex_fq_files' is provided then 'umi_regexp' should
        extract both barcodes and UMIs (i.e. it should contain both
        '<cell>' and '<umi>' elements).
    * 'dedup_umis': Deduplicate reads using UMI-tools? (default
      'FALSE')
    * 'dedup_stats': Output UMI deduplication statistics? (default
      'TRUE')
    * 'group_umis': Smmarise UMI groups both pre- and
      post-deduplication using UMI-tools? Useful for debugging
      (default 'FALSE')
    * If 'dedup_umis' is 'TRUE' but 'extract_umis' is 'FALSE' then a
      warning will be displayed, but processing will continue.
    * 'trim_5p_mismatches': Trim mismatched 5' base? (default 'TRUE')

    Statistics and figure generation input files:

    * 'asite_disp_length_file': Summary of read frame displacement
      from 5' end to A-site for each read length based on 'standard'
      yeast data from early ribosome profiling papers (tab-separated
      values file with 'read_length', 'asite_disp' columns) (optional)
    * 'codon_positions_file': Position of codons within each gene
      (RData file) (optional)
    * 'features_file': Features to correlate with ORFs (tab-separated
      values file with 'ORF', 'Length_log10', 'uATGs', 'FE_atg',
      'FE_cap', 'utr', 'utr_gc', 'polyA' columns) (optional)
    * 't_rna_file': tRNA estimates file (tab-separated values file
      with 'AA', 'Codon', 'tRNA', 'tAI', 'Microarray', 'RNA.seq'
      columns)  (optional)
    * While both 'codon_positions_file' and 't_rna_file' are optional
      either both must be specified or neither must be specified.

    Statistics and figure generation parameters:

    * 'buffer': Length of flanking region around the CDS (default 250)
    * 'count_reads': Scan input, temporary and output files and
      produce counts of reads in each FASTQ, SAM, and BAM file
      processed? (default 'TRUE')
    * 'count_threshold': Remove genes with a read count below this
      threshold, when generating statistics and figures (default 1)
    * 'dataset': Human-readable name of the dataset (default
       'dataset')
    * 'do_pos_sp_nt_freq': Calculate position-specific nucleotide
      freqeuency? (default 'TRUE')
    * 'is_riboviz_gff': Does the GFF file contain 3 elements per gene
      - UTR5, CDS, and UTR3? (default 'TRUE')
    * 'make_bedgraph': Output bedgraph data files in addition to H5
      files? (default 'TRUE')
    * 'max_read_length': Maximum read length in H5 output (default 50)
    * 'min_read_length': Minimum read length in H5 output (default 10)
    * 'primary_id': Primary gene IDs to access the data (YAL001C,
      YAL003W, etc.) (default 'Name')
    * 'rpf': Is the dataset an RPF or mRNA dataset? (default 'TRUE')
    * 'secondary_id': Secondary gene IDs to access the data (COX1,
      EFB1, etc. or 'NULL') (default 'NULL')
    * 'stop_in_cds': Are stop codons part of the CDS annotations in
      GFF? (default 'FALSE')

    General:

    * 'validate_only': Validate configuration, check that mandatory
      parameters have been provided and that input files exist, then
      exit without running the workflow? (default 'FALSE')
    * 'publish_index_tmp': Publish/copy index and temporary files to
      'dir_index' and 'dir_tmp'. If 'FALSE' then only symbolic links
      to these files in the Nextflow 'work/' directory are
      created in 'dir_index' and 'dir_tmp' (default 'FALSE')
    * 'skip_inputs': When validating configuration (see
      'validate_only' above) skip checks for existence of ribosome
      profiling data files ('fq_files', 'multiplexed_fq_files',
      'sample_sheet')? (default 'FALSE')
    * 'num_processes': Number of processes to parallelize over, used
      by specific steps in the workflow (default 1)
    * 'samsort_memory': Memory to give to 'samtools sort' (
      default '768M', 'samtools sort' built-in default,
      see http://www.htslib.org/doc/samtools-sort.html)
    """.stripIndent()
}

// Help message implementation, following
// https://github.com/nf-core/rnaseq/blob/master/main.nf (MIT License)
params.help = false
if (params.help) {
    helpMessage()
    exit 0
}

/*
Initialise and validate configuration.

Initialise optional variables to avoid "WARN: Access to undefined
parameter '<PARAM>'" errors.
*/

params.buffer = 250
params.build_indices = true
params.count_reads = true
params.count_threshold = 1
params.dataset = "dataset"
params.dedup_umis = false
params.dir_index = "index"
params.dir_out = "output"
params.dir_tmp = "tmp"
params.do_pos_sp_nt_freq = true
params.extract_umis = false
params.trim_5p_mismatches = true
params.fq_files = [:]
params.group_umis = false
params.dedup_stats = true
params.is_riboviz_gff = true
params.make_bedgraph = true
params.max_read_length = 50
params.min_read_length = 10
params.multiplex_fq_files = []
params.num_processes = 1
params.publish_index_tmp = false
params.primary_id = "Name"
params.rpf = true
params.secondary_id = "NULL"
params.stop_in_cds = false
params.samsort_memory = null
params.validate_only = false
params.skip_inputs = false

if (params.publish_index_tmp)
{
     publish_index_tmp_type = 'copy'
}
else
{
    publish_index_tmp_type = 'symlink'
}

if (params.validate_only) {
    println("Validating configuration only")
}

if (! params.containsKey('adapters')) {
    exit 1, "Undefined adapters (adapters)"
}
if (! params.containsKey('orf_index_prefix')) {
    exit 1, "Undefined ORF index prefix (orf_index_prefix)"
}
if (! params.containsKey('rrna_index_prefix')) {
    exit 1, "Undefined rRNA index prefix (rrna_index_prefix)"
}
if (params.buffer < 0) {
    exit 1, "CDS flanking region length (buffer) must be >= 0"
}
if (params.count_threshold < 0) {
    exit 1, "Read count threshold (count_threshold) must be >= 0"
}
if (params.num_processes < 1) {
    exit 1, "Number of processes (num_processes) must be >= 1"
}
if (params.min_read_length < 1) {
    exit 1, "Minimum read length in H5 output (min_read_length) must be >= 1"
}
if (params.max_read_length < 1) {
    exit 1, "Maximum read length in H5 output (max_read_length) must be >= 1"
}
if (params.max_read_length < params.min_read_length) {
    exit 1, "Maximum read length in H5 output (max_read_length) must be >= minimum read length (min_read_length)"
}
if (! params.secondary_id) {
    secondary_id = "NULL"
} else {
    secondary_id = params.secondary_id
}
if (params.dedup_umis) {
    if (! params.extract_umis) {
        println("Warning: deduplication was requested (dedup_umi: TRUE) but UMI extraction was not (extract_umis: FALSE)")
    }
}
if (params.extract_umis) {
    if (! params.containsKey('umi_regexp')) {
        exit 1, "Undefined barcode/UMI regular expression (umi_regexp) when UMI extraction is requested (extract_umis: TRUE)"
    }
}

/*
Validate input files.
*/

if (params.download) {
  to_download = []
  for (entry in params.download){
    to_download.add([entry.key,entry.value])
  }
  if (params.multiplex_fq_files){
    is_multiplexed = true
    if (! params.containsKey('sample_sheet')) {
      exit 1, "Undefined sample sheet (sample_sheet)"
    }
    sample_sheet = file("${params.dir_in}/${params.sample_sheet}")
    if ((! params.validate_only) || (! params.skip_inputs)) {
        if (! sample_sheet.exists()) {
            exit 1, "No such sample sheet (sample_sheet): ${sample_sheet}"
        }
        multiplex_sample_sheet_tsv = Channel.fromPath(sample_sheet,
                                                      checkIfExists: true)
    }
  } else{
    is_multiplexed = false
  }
  
} else {

  sample_id_fq = []
  multiplex_id_fq = []

  multiplex_sample_sheet_tsv = Channel.empty()
  is_multiplexed = false
  if (params.validate_only && params.skip_inputs) {
      println("Skipping checks for existence of of ribosome profiling input files (fq_files|multiplex_fq_files|sample_sheet)")
  }
  if (! params.containsKey('dir_in')) {
      exit 1, "Undefined input directory (dir_in)"
  } else if ((! params.fq_files) && (! params.multiplex_fq_files)) {
      exit 1, "No sample files (fq_files), multiplexed files (multiplex_fq_files)"
  } else if (params.fq_files && params.multiplex_fq_files) {
      exit 1, "Both sample files (fq_files) and multiplexed files (multiplex_fq_files) are defined - only one or the other should be defined"
  } else if (params.fq_files) {
      if ((! params.validate_only) || (! params.skip_inputs)) {
          // Filter 'params.fq_files' down to those samples that exist.
          for (entry in params.fq_files) {
              sample_fq = file("${params.dir_in}/${entry.value}")
              if (sample_fq.exists()) {
                  //sample_id_fq[entry.key] = sample_fq
                  sample_id_fq.add([entry.key,sample_fq])
              } else {
                println("No such sample file ($entry.key): $entry.value")
              }
          }
          if (! sample_id_fq) {
            exit 1, "None of the defined sample files (fq_files) exist"
          }
      }
  } else if (params.multiplex_fq_files) {
      if ((! params.validate_only) || (! params.skip_inputs)) {
          // Filter 'params.multiplex_fq_files' down to those files that exist.
          for (entry in params.multiplex_fq_files) {
              multiplex_fq = file("${params.dir_in}/${entry}")
              if (multiplex_fq.exists()) {
                  // Use file base name as key, ensuring that if file
                  // has extension '.fastq.gz' or '.fq.gz' then both
                  // extensions are removed from the name.
                  multiplex_id = multiplex_fq.baseName
                  if (multiplex_id.endsWith(".fastq")) {
                      multiplex_id = multiplex_id - '.fastq'
                  } else if (multiplex_id.endsWith(".fq")) {
                      multiplex_id = multiplex_id - '.fq'
                  }
                  multiplex_id_fq.add([multiplex_id,multiplex_fq])
              } else {
                  println("No such multiplexed file: $entry")
              }
          }
          if (! multiplex_id_fq) {
              exit 1, "None of the defined multiplexed files (multiplex_fq_files) exist"
          }
      }
      if (! params.containsKey('sample_sheet')) {
          exit 1, "Undefined sample sheet (sample_sheet)"
      }
      sample_sheet = file("${params.dir_in}/${params.sample_sheet}")
      if ((! params.validate_only) || (! params.skip_inputs)) {
          if (! sample_sheet.exists()) {
              exit 1, "No such sample sheet (sample_sheet): ${sample_sheet}"
          }
          multiplex_sample_sheet_tsv = Channel.fromPath(sample_sheet,
                                                        checkIfExists: true)
      }
      is_multiplexed = true
  } 
}


process SRA_IDS_TO_RUNINFO {
    publishDir "${params.dir_tmp}/", \
          mode: 'copy', overwrite: true
    input:
    tuple val(sample_name), val(run_id) from to_download//.collect{ sample_name, run -> [sample_name, run] }

    output:
    tuple val(sample_name), path("*.tsv") into runinfo

    when:
    params.download

    script:
    """
    ${workflow.projectDir}/riboviz/sra_ids_to_runinfo.py ${run_id} ${run_id}.runinfo.tsv
    """
  }

  process SRA_RUNINFO_TO_FTP {
      publishDir "${params.dir_tmp}", \
        mode: 'copy', overwrite: true

          
      input:
      tuple val(sample_name), path(run_id) from runinfo
      
      output:
      //Send output to both channels, but only one will run.
      path("*.tsv") into (single_fastq_info,multiplex_fastq_info) 
      val(sample_name) into (single_sample_info,multiplex_sample_info)

      script:
      """
      ${workflow.projectDir}/riboviz/sra_runinfo_to_ftp.py ${run_id.join(',')} ${run_id.toString().tokenize(".")[0]}.runinfo_ftp.tsv
      """
  }




  process SRA_FASTQ_FTP {
      publishDir "${params.dir_in}", \
        mode: 'copy', overwrite: true
      
      input:
      tuple val(single), val(fastq) from single_fastq_info.splitCsv(header:true,sep:'\t').map{ row -> [row.single_end,row.fastq_1]}
      val(sample_name) from single_sample_info

      output:
      
      tuple val(sample_name), path("*fastq.gz") optional true into sample_id_fq

      when:
      (! is_multiplexed)
      script:
      // Based on call to hisat2, seems like we are only concerned with single-end reads. For now, don't worry about paired-end.
      if (single.toBoolean()){
        """
        curl -L ${fastq[0]} -O
        """
      } else{
        """
        echo 'WARNING: Provided SRA dataset is paired-end read data. Riboviz currently only handles single-end read data. These reads will not be downloaded.'
        """
      }

      
  }

/*TODO: @acope3 I don't know if we need this separate process for multiplexed data. I think we can just have SRA_FASTQ_FTP output 
to both sample_id_fq and multiplex_id_fq channels, with the appropriate cutAdapt process being called based on value of is_multiplexed.
*/
process SRA_FASTQ_FTP_Multiplexed {
      publishDir "${params.dir_in}", \
        mode: 'copy', overwrite: true
      
      input:
      tuple val(single), val(fastq) from multiplex_fastq_info.splitCsv(header:true,sep:'\t').map{ row -> [row.single_end,row.fastq_1]}
      val(sample_name) from multiplex_sample_info

      output:
      tuple val(sample_name), path("*fastq.gz") into multiplex_id_fq

      when:
      (is_multiplexed)
      script:
      // Based on call to hisat2, seems like we are only concerned with single-end reads. For now, don't worry about paired-end.
      if (single.toBoolean()){
        """
        curl -L ${fastq[0]} -O
        """
      } else{
        """
        echo 'WARNING: Provided SRA dataset is paired-end read data. Riboviz currently only handles single-end read data. These reads will not be downloaded.'
        """
      }
      
  }

// Create YAML fragment including 'params.fq_files' and
// 'params.multiplex_fq_files' to serve as a configuration file for
// riboviz.tools.count_reads. There is no way to get the location
// of the YAML configuration file itself (i.e. the value of
// '-param-file') from within Nextflow so this is a workaround.
Map ribosome_fqs = [:]
ribosome_fqs.fq_files = params.fq_files
ribosome_fqs.multiplex_fq_files = params.multiplex_fq_files
ribosome_fqs_yaml = new Yaml().dump(ribosome_fqs)

// Non-sample-specific input files.
if (! params.build_indices) {
    rrna_index_prefix = file("${params.dir_index}/${params.rrna_index_prefix}.*.ht2")
    if (! rrna_index_prefix) {
        exit 1, "No such rRNA index files (rrna_index_prefix): ${params.dir_index}/${params.rrna_index_prefix}.*.ht2"
    }
    pre_built_rrna_index_ht2 = Channel
        .fromPath(rrna_index_prefix, checkIfExists: true)
        .collect()
    orf_index_prefix = file("${params.dir_index}/${params.orf_index_prefix}.*.ht2")
    if (! orf_index_prefix) {
        exit 1, "No such ORF index files (orf_index_prefix): ${params.dir_index}/${params.orf_index_prefix}.*.ht2"
    }
    pre_built_orf_index_ht2 = Channel
        .fromPath("${params.dir_index}/${params.orf_index_prefix}.*.ht2",
                  checkIfExists: true)
        .collect()
} else {
    pre_built_rrna_index_ht2 = Channel.empty()
    pre_built_orf_index_ht2 = Channel.empty()
}

if (! params.containsKey('rrna_fasta_file')) {
    exit 1, "Undefined rRNA FASTA file (rrna_fasta_file)"
}
rrna_fasta_file = file(params.rrna_fasta_file)
if (! rrna_fasta_file.exists()) {
    exit 1, "No such file rRNA FASTA file (rrna_fasta_file): ${rrna_fasta_file}"
}
rrna_fasta = Channel.fromPath(rrna_fasta_file, checkIfExists: true)
if (! params.containsKey('orf_fasta_file')) {
    exit 1, "Undefined ORF FASTA file (orf_fasta_file)"
}
orf_fasta_file = file(params.orf_fasta_file)
if (! orf_fasta_file.exists()) {
    exit 1, "No such ORF FASTA file (orf_fasta_file): ${orf_fasta_file}"
}
orf_fasta = Channel.fromPath(orf_fasta_file, checkIfExists: true)
if (! params.containsKey('orf_gff_file')) {
    exit 1, "Undefined ORF GFF file (orf_gff_file)"
}
orf_gff_file = file(params.orf_gff_file)
if (! orf_gff_file.exists()) {
    exit 1, "No such ORF GFF file (orf_gff_file): ${orf_gff_file}"
}
orf_gff = Channel.fromPath(orf_gff_file, checkIfExists: true)

// Optional inputs for generate_stats_figs.R.
// If an optional file is not provided then a 'Missing_<PARAM>' file
// (for example 'Missing_features_file') is created within the 'work/'
// directories for the generateStatsFigs process. This symbolically
// links to a non-existent 'Missing_<PARAM>' file in the users current
// directory. This is not an issue since the files will not be passed
// onto generate_stats_figs.R and no attempt is made to use them. They
// are a side-effect of using the Nextflow pattern for optional
// inputs.
// Optional inputs implementation follows pattern
// https://github.com/nextflow-io/patterns/blob/master/optional-input.nf.
if (params.containsKey('t_rna_file') && params.t_rna_file) {
    t_rna_file = file(params.t_rna_file)
    if (! t_rna_file.exists()) {
        exit 1, "No such tRNA estimates file (t_rna_file): ${t_rna_file}"
    }
    t_rna_tsv = Channel.fromPath(t_rna_file, checkIfExists: true)
    is_t_rna_file = true
} else {
    t_rna_tsv = file("Missing_t_rna_file")
    is_t_rna_file = false
}
if (params.containsKey('codon_positions_file')
    && params.codon_positions_file) {
    codon_positions_file = file(params.codon_positions_file)
    if (! codon_positions_file.exists()) {
        exit 1, "No such codon positions file (codon_positions_file): ${codon_positions_file}"
    }
    codon_positions_rdata = Channel.fromPath(codon_positions_file,
                                             checkIfExists: true)
    is_codon_positions_file = true
} else {
    codon_positions_rdata = file("Missing_codon_positions_file")
    is_codon_positions_file = false
}
if (is_t_rna_file && is_codon_positions_file) {
    is_t_rna_and_codon_positions_file = true
} else if ((! is_t_rna_file) && (! is_codon_positions_file)) {
    is_t_rna_and_codon_positions_file = false
} else {
    exit 1, "Either both tRNA estimates (t_rna_file) and codon positions (codon_positions_file) must be defined or neither must be defined"
}
if (params.containsKey('features_file') && params.features_file) {
    features_file = file(params.features_file)
    if (! features_file.exists()) {
        exit 1, "No such features file (features_file): ${features_file}"
    }
    features_tsv = Channel.fromPath(features_file, checkIfExists: true)
    is_features_file = true
} else {
    features_tsv = file("Missing_features_file")
    is_features_file = false
}
if (params.containsKey('asite_disp_length_file')
    && params.asite_disp_length_file) {
    asite_disp_length_file = file(params.asite_disp_length_file)
    if (! asite_disp_length_file.exists()) {
        exit 1, "No such A-site displacement file (asite_disp_length_file): ${asite_disp_length_file}"
    }
    asite_disp_length_txt = Channel.fromPath(asite_disp_length_file,
                                             checkIfExists: true)
    is_asite_disp_length_file = true
} else {
    asite_disp_length_txt = file("Missing_aside_disp_length_file")
    is_asite_disp_length_file = false
}

if (params.validate_only) {
    exit 0, "Validated configuration"
}

/*
Indexing.
*/

// Split channels for use in multiple downstream processes.
orf_fasta.into { build_indices_orf_fasta; generate_stats_figs_orf_fasta }
orf_gff.into { bam_to_h5_orf_gff; generate_stats_figs_orf_gff }

process buildIndicesrRNA {
    tag "${params.rrna_index_prefix}"
    publishDir "${params.dir_index}", mode: publish_index_tmp_type, overwrite: true
    input:
        file rrna_fasta from rrna_fasta
    output:
        file "${params.rrna_index_prefix}.*.ht2" into built_rrna_index_ht2
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${rrna_fasta} ${params.rrna_index_prefix}
        """
}

process buildIndicesORF {
    tag "${params.orf_index_prefix}"
    publishDir "${params.dir_index}", mode: publish_index_tmp_type, overwrite: true
    input:
        file orf_fasta from build_indices_orf_fasta
    output:
        file "${params.orf_index_prefix}.*.ht2" into built_orf_index_ht2
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${orf_fasta} ${params.orf_index_prefix}
        """
}

// Combine "[pre_]built_rrna|orf_index_ht2" channels for downstream
// processing. Due to foregoing behaviour conditional on
// 'params.build_indices' one one of 'pre_built_rrna|orf_index_ht2' or
// 'built_rrna|orf_index_ht2' will have content.
rrna_index_ht2 = pre_built_rrna_index_ht2.mix(built_rrna_index_ht2)
orf_index_ht2 = pre_built_orf_index_ht2.mix(built_orf_index_ht2)

/*
Sample file (fq_files)-specific processes.
*/

process cutAdapters {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(sample_id), file(sample_fq) \
            from sample_id_fq//.collect{ id, file -> [id, file] }
    output:
        tuple val(sample_id), file("trim.fq") into cut_fq
    when:
        (! is_multiplexed)
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o trim.fq ${sample_fq} -j 0
        """
}

// Route 'cut_fq' channel outputs depending on whether UMIs are to be
// extracted or not.
cut_fq.branch {
    umi_fq: params.extract_umis
    non_umi_fq: ! params.extract_umis
}
.set { cut_fq_branch }

process extractUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(sample_id), file(sample_fq) \
            from cut_fq_branch.umi_fq
    output:
        tuple val(sample_id), file("extract_trim.fq") \
            into umi_extract_fq
    when:
        params.extract_umis && (! is_multiplexed)
    shell:
        """
        umi_tools extract -I ${sample_fq} \
            --bc-pattern="${params.umi_regexp}" \
            --extract-method=regex -S extract_trim.fq
        """
}

/*
Multiplexed files (multiplex_fq_files)-specific processes.
*/

process cutAdaptersMultiplex {
    tag "${multiplex_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}", mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(multiplex_id), file(multiplex_fq) \
            from multiplex_id_fq//.collect{ id, file -> [id, file] }
    output:
        tuple val(multiplex_id), file("${multiplex_id}_trim.fq") \
            into cut_multiplex_fq
    when:
        is_multiplexed
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o ${multiplex_id}_trim.fq ${multiplex_fq} -j 0
        """
}

// Route 'cut_multiplex_fq' channel outputs depending on whether UMIs
// are to be extracted or not.
cut_multiplex_fq.branch {
    umi_fq: params.extract_umis
    non_umi_fq: ! params.extract_umis
}
.set { cut_multiplex_fq_branch }

process extractUmisMultiplex {
    tag "${multiplex_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}", mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(multiplex_id), file(multiplex_fq) \
            from cut_multiplex_fq_branch.umi_fq
    output:
        tuple val(multiplex_id), file("${multiplex_id}_extract_trim.fq") \
            into umi_extract_multiplex_fq
    when:
        params.extract_umis && is_multiplexed
    shell:
        """
        umi_tools extract -I ${multiplex_fq} \
            --bc-pattern="${params.umi_regexp}" \
            --extract-method=regex -S ${multiplex_id}_extract_trim.fq
        """
}

// Combine channels for downstream processing. By definition of
// 'cut_multiplex_fq.branch' only one of the input channels will have
// content.
trimmed_multiplex_fq = cut_multiplex_fq_branch.non_umi_fq
    .mix(umi_extract_multiplex_fq)

// Split channel for use in multiple downstream processes.
multiplex_sample_sheet_tsv.into {
    report_multiplex_sample_sheet_tsv; deplex_multiplex_sample_sheet_tsv
}

process demultiplex {
    tag "${multiplex_id}"
    publishDir "${params.dir_tmp}/${multiplex_id}_deplex", \
        mode: publish_index_tmp_type, overwrite: true
    errorStrategy 'ignore'
    input:
        // Use '.toString' to prevent changing hashes of
        // 'workflow.projectDir' triggering reexecution of this
        // process if 'nextflow run' is run with '-resume'.
        env PYTHONPATH from workflow.projectDir.toString()
        tuple val(multiplex_id), file(multiplex_fq) from trimmed_multiplex_fq
        each file(sample_sheet_tsv) from deplex_multiplex_sample_sheet_tsv
    output:
        tuple val(multiplex_id), file("num_reads.tsv") \
                into demultiplex_num_reads_tsv
        file("*.f*") into demultiplex_fq
    shell:
        """
        python -m riboviz.tools.demultiplex_fastq \
            -1 ${multiplex_fq} -s ${sample_sheet_tsv} -o . -m 2
        """
}

// 'demultiplex_fq' outputs a single list with all the output
// files. Extract sample IDs from file basenames, filter out
// 'Unassigned' and output tuples of sample IDs and file names as
// separate items onto a new channel.
demultiplex_fq
    .flatten()
    // Use file basename as sample ID.
    .map { [it.baseName, it] }
    // If file was '.fastq|fq.gz' then basename will include
    // '.fq|fastq' so strip that off too.
    .map { n, f -> [n.endsWith(".fq") ? n - ".fq" : n, f] }
    .map { n, f -> [n.endsWith(".fastq") ? n - ".fastq" : n, f] }
    .filter { n, f -> n != "Unassigned" }
    .into { report_demultiplex_samples_fq; demultiplex_samples_fq }

if (is_multiplexed) {

    demultiplex_sample_ids = report_demultiplex_samples_fq
        .map { n, f -> n }
        .toList() // [] if none
        .view { "Demultiplexed samples: ${it}"}
        // Wrap list in list, so 'merge' below doesn't append lists
        .map { it -> [it] }

    multiplex_sample_sheet_ids = report_multiplex_sample_sheet_tsv
        // Extract original sample IDs from sample sheet
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.SampleID } // No output if no 'SampleID' column
        .toList() // [] if no 'SampleID' column
        // Wrap list in list, so 'merge' below doesn't append lists
        .map { it -> [it] }

    multiplex_sample_sheet_ids
        .merge(demultiplex_sample_ids)
        .map { a, b -> a - b}
        .view { "Non-demultiplexed samples: ${it}" }
}

/*
Sample-specific processes.

Common to both sample files (fq_files) and demultiplexed files.
*/

// Combine channels for downstream processing. By definition of
// upstream conditions and processes, only one of the channels
// will have content.
trimmed_fq = cut_fq_branch.non_umi_fq
    .mix(umi_extract_fq)
    .mix(demultiplex_samples_fq)

process hisat2rRNA {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sample_fq) from trimmed_fq
        each file(rrna_index_ht2) from rrna_index_ht2
    output:
        tuple val(sample_id), file("nonrRNA.fq") into non_rrna_fq
        tuple val(sample_id), file("rRNA_map.sam") into rrna_map_sam
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -N 1 -k 1 \
            --un nonrRNA.fq -x ${params.rrna_index_prefix} \
            -S rRNA_map.sam -U ${sample_fq}
        """
}

process hisat2ORF {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sample_fq) from non_rrna_fq
        each file(orf_index_ht2) from orf_index_ht2
    output:
        tuple val(sample_id), file("unaligned.fq") into unaligned_fq
        tuple val(sample_id), file("orf_map.sam") into trim_5p_mismatches
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -k 2 \
            --no-spliced-alignment --rna-strandness F --no-unal \
            --un unaligned.fq -x ${params.orf_index_prefix} \
            -S orf_map.sam -U ${sample_fq}
        """
}

// Route 'trim_5p_branch' channel outputs depending on whether mismatched
// 5' base are to be trimmed or not
trim_5p_mismatches.branch {
    trim_5p_fq: params.trim_5p_mismatches
    non_trim_5p_fq: ! params.trim_5p_mismatches
}
.set { trim_5p_branch }

process trim5pMismatches {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    errorStrategy 'ignore'
    input:
        // Use '.toString' to prevent changing hashes of
        // 'workflow.projectDir' triggering reexecution of this
        // process if 'nextflow run' is run with '-resume'.
        env PYTHONPATH from workflow.projectDir.toString()
        tuple val(sample_id), file(sample_sam) from trim_5p_branch.trim_5p_fq
    output:
        tuple val(sample_id), file("orf_map_clean.sam") \
            into trim_orf_map_sam
        tuple val(sample_id), file("trim_5p_mismatch.tsv") \
            into trim_summary_tsv
    shell:
        """
        python -m riboviz.tools.trim_5p_mismatch -m 2 \
            -i ${sample_sam} -o orf_map_clean.sam -s trim_5p_mismatch.tsv
        """
}

// Combine channels for downstream processing. By definition of
// upstream conditions and processes, only one of the channels
// will have content.
trimmed_5p_fq = trim_5p_branch.non_trim_5p_fq
    .mix(trim_orf_map_sam)

process samViewSort {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sample_sam) from trimmed_5p_fq
    output:
        tuple val(sample_id), file("orf_map_clean.bam"), \
            file("orf_map_clean.bam.bai") into orf_map_bam
    shell:
        memory = params.samsort_memory != null ? "-m ${params.samsort_memory}" : ""
        """
        samtools --version
        samtools view -b ${sample_sam} | samtools sort ${memory} \
            -@ ${params.num_processes} -O bam -o orf_map_clean.bam -
        samtools index orf_map_clean.bam
        """
}

// Route "orf_map_bam" channel outputs depending on whether UMIs are
// to be deduplicated or not.
orf_map_bam.branch {
    dedup_bam: params.dedup_umis
    non_dedup_bam: ! params.dedup_umis
}
.set { orf_map_bam_branch }

// Split channel for use in multiple downstream processes.
orf_map_bam_branch.dedup_bam.into {
    pre_dedup_group_bam; pre_dedup_bam
}

process groupUmisPreDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(sample_id), file(sample_bam), file(sample_bam_bai) \
            from pre_dedup_group_bam
    output:
        tuple val(sample_id), file("pre_dedup_groups.tsv") \
            into pre_dedup_group_tsv
    when:
        params.dedup_umis && params.group_umis
    shell:
        """
        umi_tools group -I ${sample_bam} --group-out pre_dedup_groups.tsv
        """
}

process dedupUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(sample_id), file(sample_bam), file(sample_bam_bai) \
            from pre_dedup_bam
    output:
        tuple val(sample_id), file("dedup.bam"), \
            file("dedup.bam.bai") into dedup_bam
        tuple val(sample_id), file("dedup_stats*.tsv") \
            optional (! params.dedup_stats) \
            into dedup_stats_tsv
    when:
        params.dedup_umis
    shell:
        output_stats_flag = params.dedup_stats \
	    ? "--output-stats=dedup_stats" : ''
        """
        umi_tools dedup -I ${sample_bam} -S dedup.bam ${output_stats_flag}
        samtools --version
        samtools index dedup.bam
        """
}

// Split channel for use in multiple downstream processes.
dedup_bam.into { post_dedup_group_bam; post_dedup_bam }

process groupUmisPostDedup {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${params.dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(sample_id), file(sample_bam), file(sample_bam_bai) \
            from post_dedup_group_bam
    output:
        tuple val(sample_id), file("post_dedup_groups.tsv") \
            into post_dedup_group_tsv
    when:
        params.dedup_umis && params.group_umis
    shell:
        """
        umi_tools group -I ${sample_bam} --group-out post_dedup_groups.tsv
        """
}

// Combine channels for downstream processing. By definition of
// 'orf_map_bam_branch' only one of the input channels will have
// content.
pre_output_bam = orf_map_bam_branch.non_dedup_bam.mix(post_dedup_bam)

process outputBams {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sample_bam), file(sample_bam_bai) \
            from pre_output_bam
    output:
        tuple val(sample_id), file("${sample_id}.bam"), \
            file("${sample_id}.bam.bai") into output_bam
    shell:
        """
        cp ${sample_bam} ${sample_id}.bam
        cp ${sample_bam_bai} ${sample_id}.bam.bai
        """
}

// Split channel for use in multiple downstream processes.
output_bam.into { bedgraph_bam; bam_to_h5_bam }

process makeBedgraphs {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sample_bam), file(sample_bam_bai) \
            from bedgraph_bam
    output:
        tuple val(sample_id), file("plus.bedgraph"), \
            file("minus.bedgraph") into bedgraph
    when:
        params.make_bedgraph
    shell:
        """
        bedtools --version
        bedtools genomecov -ibam ${sample_bam} -trackline -bga -5 \
            -strand + > plus.bedgraph
        bedtools genomecov -ibam ${sample_bam} -trackline -bga -5 \
            -strand - > minus.bedgraph
        """
}

process bamToH5 {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sample_bam), \
            file(sample_bam_bai) from bam_to_h5_bam
        each file(orf_gff) from bam_to_h5_orf_gff
    output:
        tuple val(sample_id), file("${sample_id}.h5") into h5s
    shell:
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/bam_to_h5.R \
           --num-processes=${params.num_processes} \
           --min-read-length=${params.min_read_length} \
           --max-read-length=${params.max_read_length} \
           --buffer=${params.buffer} \
           --primary-id=${params.primary_id} \
           --secondary-id=${secondary_id} \
           --dataset=${params.dataset} \
           --bam-file=${sample_bam} \
           --hd-file=${sample_id}.h5 \
           --orf-gff-file=${orf_gff} \
           --is-riboviz-gff=${params.is_riboviz_gff} \
           --stop-in-cds=${params.stop_in_cds}
        """
}

// Optional inputs implementation follows pattern
// https://github.com/nextflow-io/patterns/blob/master/optional-input.nf.
process generateStatsFigs {
    tag "${sample_id}"
    publishDir "${params.dir_out}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(sample_h5) from h5s
        each file(orf_fasta) from generate_stats_figs_orf_fasta
        each file(orf_gff) from generate_stats_figs_orf_gff
        each file(t_rna_tsv) from t_rna_tsv
        each file(codon_positions_rdata) from codon_positions_rdata
        each file(features_tsv) from features_tsv
        each file(asite_disp_length_txt) from asite_disp_length_txt
    output:
        val sample_id into finished_sample_id
        tuple val(sample_id), file("tpms.tsv") into tpms_tsv
        tuple val(sample_id), file("3nt_periodicity.pdf") \
            into nt3_periodicity_pdf
        tuple val(sample_id), file("3nt_periodicity.tsv") \
            into nt3_periodicity_tsv
        tuple val(sample_id), file("pos_sp_nt_freq.tsv") \
	    optional (! params.do_pos_sp_nt_freq) \
	    into pos_sp_nt_freq_tsv
        tuple val(sample_id), file("pos_sp_rpf_norm_reads.pdf") \
            into pos_sp_rpf_norm_reads_pdf
        tuple val(sample_id), file("pos_sp_rpf_norm_reads.tsv") \
            into pos_sp_rpf_norm_reads_tsv
        tuple val(sample_id), file("read_lengths.pdf") \
            into read_lengths_pdf
        tuple val(sample_id), file("read_lengths.tsv") \
            into read_lengths_tsv
        tuple val(sample_id), file("startcodon_ribogridbar.pdf") \
            into start_codon_ribogridbar_pdf
        tuple val(sample_id), file("startcodon_ribogrid.pdf") \
            into start_codon_ribogrid_pdf
        tuple val(sample_id), file("codon_ribodens.pdf") \
            optional (! is_t_rna_and_codon_positions_file) \
            into codon_ribodens_pdf
        tuple val(sample_id), file("codon_ribodens.tsv") \
            optional (! is_t_rna_and_codon_positions_file) \
            into codon_ribodens_tsv
        tuple val(sample_id), file("features.pdf") \
            optional (! is_features_file) into features_pdf
        tuple val(sample_id), file("3ntframe_bygene.tsv") \
            optional (! is_asite_disp_length_file) \
            into nt3frame_bygene_tsv
        tuple val(sample_id), file("3ntframe_propbygene.pdf") \
            optional (! is_asite_disp_length_file) \
            into nt3frame_propbygene_pdf
    shell:
        t_rna_flag = is_t_rna_and_codon_positions_file \
            ? "--t-rna-file=${t_rna_tsv}" : ''
        codon_positions_flag = is_t_rna_and_codon_positions_file \
            ? "--codon-positions-file=${codon_positions_rdata}" : ''
        features_flag = is_features_file \
            ? "--features-file=${features_tsv}" : ''
        asite_disp_length_flag = is_asite_disp_length_file \
            ? "--asite-disp-length-file=${asite_disp_length_txt}" : ''
        count_threshold_flag = params.containsKey('count_threshold') \
            ? "--count-threshold=${params['count_threshold']}": ''
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/generate_stats_figs.R \
           --num-processes=${params.num_processes} \
           --min-read-length=${params.min_read_length} \
           --max-read-length=${params.max_read_length} \
           --buffer=${params.buffer} \
           --primary-id=${params.primary_id} \
           --dataset=${params.dataset} \
           --hd-file=${sample_h5} \
           --orf-fasta-file=${orf_fasta} \
           --rpf=${params.rpf} \
           --output-dir=. \
           --do-pos-sp-nt-freq=${params.do_pos_sp_nt_freq} \
           ${t_rna_flag} \
           ${codon_positions_flag} \
           ${features_flag} \
           --orf-gff-file=${orf_gff} \
           ${asite_disp_length_flag} \
           ${count_threshold_flag}
        """
}

finished_sample_id
    .ifEmpty { exit 1, "No sample was processed successfully" }
    .view { "Finished processing sample: ${it}" }

// Prefix sample-specific TPMs files, tpms.tsv, with sample ID so all
// sample-specific TPMs files can be staged into the same directory
// for running collateTpms.
process renameTpms {
    tag "${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), file(tpms_tsv) from tpms_tsv
    output:
        val(sample_id) into tpms_sample_id
        file "${sample_id}_tpms.tsv" into tpms_sample_tsv
    shell:
        """
        cp ${tpms_tsv} ${sample_id}_tpms.tsv
        """
}

process collateTpms {
    tag "${sample_ids.join(', ')}"
    publishDir "${params.dir_out}", mode: 'copy', overwrite: true
    input:
        val sample_ids from tpms_sample_id.collect()
        file tpms_tsvs from tpms_sample_tsv.collect()
    output:
        file "TPMs_collated.tsv" into collate_tpms_tsv
        val sample_ids into collate_tpms_sample_ids
    shell:
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/collate_tpms.R \
            --sample-subdirs=False \
            --output-dir=. \
            --tpms-file=TPMs_collated.tsv \
            ${sample_ids.join(' ')}
        """
}

process countReads {
    publishDir "${params.dir_out}", mode: 'copy', overwrite: true
    input:
        // Use '.toString' to prevent changing hashes of
        // 'workflow.projectDir' triggering reexecution of this
        // process if 'nextflow run' is run with '-resume'.
        env PYTHONPATH from workflow.projectDir.toString()
        val ribosome_fqs_yaml from ribosome_fqs_yaml
        // Force dependency on output of 'collateTpms' so this process
        // is only run when all other processing has completed.
        val samples_ids from collate_tpms_sample_ids
    output:
        file "read_counts.tsv" into read_counts_tsv
    when:
        params.count_reads
    shell:
        // 'workflow.projectDir' is directory into which outputs
        // have been published.
        // TODO: It would be preferable to:
        // 1. Stage these directories into this process's
        //    work directory. If this is possible, and how to do it,
        //    if so, has not been determined.
        // 2. Stage outputs from the processes for which reads
        //    are to be counted into this process's work directory.
        //    This would require a new implementation of
        //    riboviz.tools.count_reads.
        """
        echo "${ribosome_fqs_yaml}" > ribosome_fqs.yaml
        python -m riboviz.tools.count_reads \
           -c ribosome_fqs.yaml \
           -i ${file(params.dir_in).toAbsolutePath()} \
           -t ${file(params.dir_tmp).toAbsolutePath()} \
           -o ${file(params.dir_out).toAbsolutePath()} \
           -r read_counts.tsv
        """
}

workflow.onComplete {
    println "Workflow finished! (${workflow.success ? 'OK' : 'failed'})"
}
