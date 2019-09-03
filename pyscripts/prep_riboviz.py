#!/usr/bin/env python
"""
RiboViz workflow.

Usage:

    python pyscripts/prepRiboViz.py \
        <PYTHON_SCRIPTS_DIRECTORY>\
        <R_SCRIPTS_DIRECTORY>\
        <DATA_DIRECTORY>\
        <YAML_CONFIG_FILE>

Example:

    python pyscripts/prepRiboviz.py pyscripts/ rscripts/ data/\
        vignette/vignette_config.yaml

Prepare ribosome profiling data for RiboViz or other analysis:

* - Read configuration info from YAML file.
* - Build hisat2 indices if requested (config["build_indices"] ==
*   True), in index directory (config["dir_index"]).
* - Process all fastq.gz files (config["dir_in"]).
* - Cut out sequencing library adapters ("CTGTAGGCACC" or
*   config["adapters").
* - Remove rRNA or other contaminating reads by hisat2 alignment to
*   rRNA index file (config["rRNA_index"])
* - Align remaining reads to ORFs or other hisat2 index file
*   (config["orf_index"]).
* - Trim 5' mismatches from reads and remove reads with more than 2
*   mismatches.
* - Parallelize over many processes (config["nprocesses"]):
*   - This value is used to configure hisat2, samtools sort,
*     bam_to_h5.R and generate_stats_figs.R.
*   - For cutadapt and Python 3, the number of available processors
*     on the host will be used.
*   - For cutadapt and Python 2, its default of 1 processor will be
*     used as cutadapt cannot run in parallel under Python 2.
* - Make length-sensitive alignments in compressed h5 format by
*   running "bam_to_h5.R".
* - Generate summary statistics, and analyses and QC plots for both
*   RPF and mRNA datasets, by running "generate_stats_figs.R".
* - Put all intermediate files into a temporary directory
*   (config["dir_tmp"]).
* - When finished, put useful output files into output directory
*   (config["dir_out"]).
* Optionally export bedgraph files for plus and minus strands, if
*   requested (config["make_bedgraph"] == True).
*  - bamfiles (config["dir_out"]/*.bam) are directly usable in genome
*  - browsers, etc.

Exit codes are as follows:

* EXIT_OK (0): Processing successfully completed.
* EXIT_CONFIG_ERROR (1): Errors occurred loading configuration.
* EXIT_INDEX_ERROR (2): Error occurred during indexing.
* EXIT_NO_SAMPLES_ERROR (3): No samples were provided.
* EXIT_SAMPLES_ERROR (4): No sample was processed successfully.
* EXIT_COLLATION_ERROR (5): Error occurred during TPMs collation.
"""

import os
import os.path
import subprocess
import sys
import traceback
import yaml


EXIT_OK = 0
""" Processing successfully completed. """
EXIT_CONFIG_ERROR = 1
""" Errors occurred loading configuration. """
EXIT_INDEX_ERROR = 2
""" Error occurred during indexing. """
EXIT_NO_SAMPLES_ERROR = 3
""" No samples were provided. """
EXIT_SAMPLES_ERROR = 4
""" No sample was processed successfully. """
EXIT_COLLATION_ERROR = 5
""" Error occurred during TPMs collation. """


def list_to_str(lst):
    """
    Convert list to space-delimited string.

    :param lst: list
    :type lst: list
    :return: list as string
    :rtype: str or unicode
    """
    return ' '.join(map(str, lst))


def run_command(cmd):
    """
    Helper function to run shell command.

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :raise FileNotFoundError: if the command being run cannot be found
    :raise AssertionError: if the command returns a non-zero exit code
    """
    print(("Running: " + list_to_str(cmd)))
    exit_code = subprocess.call(cmd)
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)


def run_redirect_command(cmd, out_file):
    """
    Helper function to run shell command and redirect output to a file.

    Use pattern suggested by:
    https://www.saltycrane.com/blog/2008/09/how-get-stdout-and-stderr-using-python-subprocess-module/

    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param out_file: Output file
    :type out_file: str or unicode
    :raise FileNotFoundError: if the command being run cannot be found
    :raise AssertionError: if the command returns a non-zero exit code
    """
    print(("Running: " + list_to_str(cmd)))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdout, _ = p.communicate()
    exit_code = p.returncode
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)
    with open(out_file, "wb") as f:
        f.write(stdout)


def run_pipe_command(cmd1, cmd2):
    """
    Helper function to run shell command and pipe output into another.

    Use pattern suggested by:
    https://docs.python.org/2/library/subprocess.html#replacing-shell-pipeline
    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :param cmd: Commnand to run and its arguments
    :type cmd: list(str or unicode)
    :raise FileNotFoundError: if the commands being run cannot be found
    :raise AssertionError: if the commands returns a non-zero exit code
    """
    print(("Running: " + list_to_str(cmd1) + " | " + list_to_str(cmd2)))
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    process2 = subprocess.Popen(cmd2, stdin=process1.stdout)
    process1.stdout.close()
    output, _= process2.communicate()
    exit_code = process2.returncode
    assert exit_code == 0, ("%s % %s failed with exit code %d"
                            % (cmd1, cmd2, exit_code))


def build_indices(fasta, ht_prefix):
    """
    Build indices for alignment via invocation of hisat2-build.
    Index files have name <ht_prefix>.<N>.ht2.

    :param fasta: FASTA file to be indexed
    :type fasta: str or unicode
    :param ht_prefix: Prefix of HT2 index files
    :type ht_prefix: str or unicode
    :raise FileNotFoundError: if hisat2-build cannot be found
    :raise AssertionError: if hisat2-build returns non-zero exit
    code
    """
    cmd = ["hisat2-build", fasta, ht_prefix]
    print(("Running: " + list_to_str(cmd)))
    exit_code = subprocess.call(cmd)
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)


def process_sample(sample,
                   fastq,
                   r_rna_index,
                   orf_index,
                   config,
                   py_scripts,
                   r_scripts,
                   data_dir,
                   tmp_dir,
                   out_dir):
    """
    Process a single FASTQ sample file.

    :param sample: Sample name
    :type sample: str or unicode
    :param fastq: Sample FASTQ file
    :type fastq: str or unicode
    :param r_rna_index: Prefix of rRNA HT2 index files
    :type r_rna_index: str or unicode
    :param orf_index: Prefix of ORF HT2 index files
    :type orf_index: str or unicode
    :param config: RiboViz configuration
    :type config: dict
    :param py_scripts: Directory with RiboViz Python scripts
    :type py_scripts: str or unicode
    :param r_scripts:  Directory with RiboViz R scripts
    :type r_scripts: str or unicode
    :param data_dir: Directory with data
    :type data_dir: str or unicode
    :param tmp_dir Temporary files directory
    :type tmp_dir: str or unicode
    :param out_dir Output files directory
    :type out_dir: str or unicode
    :raise FileNotFoundError: if fastq or a third-party tool cannot be
    found
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    :raise KeyError: if config is missing required configuration
    """
    print(("Processing sample " + sample))
    if not os.path.exists(fastq):
        print(("File " + fastq + " not found"))
        raise FileNotFoundError(fastq)
    print(("Processing file " + fastq))

    if "nprocesses" not in config:
        nprocesses = 1
    else:
        nprocesses = config["nprocesses"]
    
    # Cut illumina adapters.

    # Trimmed reads.
    trim_fq = os.path.join(tmp_dir, sample + "_trim.fq")
    cmd = ["cutadapt", "--trim-n", "-O", "1", "-m", "5",
           "-a", config["adapters"], "-o", trim_fq, fastq]
    py_major = sys.version_info.major
    if py_major == 3:
        # cutadapt and Python 3 allows all available processors to
        # be requested.
        cmd += ["-j", str(0)]
    run_command(cmd)

    # Map reads to rRNA.

    # Trimmed non-rRNA reads.
    non_r_rna_trim_fq = os.path.join(tmp_dir, sample + "_nonrRNA.fq")
    # rRNA-mapped reads.
    r_rna_map_sam = os.path.join(tmp_dir, sample + "_rRNA_map.sam")
    cmd = ["hisat2", "-p", str(nprocesses), "-N", "1",
           "--un", non_r_rna_trim_fq, "-x", r_rna_index,
           "-S", r_rna_map_sam, "-U", trim_fq]
    run_command(cmd)

    # Map to ORFs with (mostly) default settings, up to 2 alignments.

    # ORF-mapped reads.
    orf_map_sam = os.path.join(tmp_dir, sample + "_orf_map.sam")
    # Unaligned reads.
    unaligned_fq = os.path.join(tmp_dir, sample + "_unaligned.fq")
    cmd = ["hisat2", "-p", str(nprocesses), "-k", "2",
           "--no-spliced-alignment", "--rna-strandness",
           "F", "--no-unal", "--un", unaligned_fq,
           "-x", orf_index, "-S", orf_map_sam,
           "-U", non_r_rna_trim_fq]
    run_command(cmd)

    # Trim 5' mismatched nt and remove reads with >1 mismatch.

    # ORF-mapped reads.
    orf_map_sam_clean = os.path.join(tmp_dir,
                                     sample +
                                     "_orf_map_clean.sam")
    cmd = ["python", os.path.join(py_scripts, "trim_5p_mismatch.py"),
           "-mm", "2", "-in", orf_map_sam,
           "-out", orf_map_sam_clean]
    run_command(cmd)

    # Convert SAM (text) output to BAM (compressed binary) and sort on
    # genome.

    sample_out_prefix = os.path.join(out_dir, sample)
    sample_out_bam = sample_out_prefix + ".bam"
    cmd_view = ["samtools", "view", "-b", orf_map_sam_clean]
    cmd_sort = ["samtools", "sort", "-@", str(nprocesses),
                "-O", "bam", "-o", sample_out_bam, "-"]
    run_pipe_command(cmd_view, cmd_sort)

    # Index BAM file.
    cmd = ["samtools", "index", sample_out_bam]
    run_command(cmd)

    if config["make_bedgraph"]:
        # Record transcriptome coverage as a bedgraph.
        # Calculate transcriptome coverage for plus strand.
        cmd = ["bedtools", "genomecov", "-ibam", sample_out_bam,
               "-trackline", "-bga", "-5", "-strand", "+"]
        run_redirect_command(cmd, sample_out_prefix + "_plus.bedgraph")
        # Calculate transcriptome coverage for minus strand.
        cmd = ["bedtools", "genomecov", "-ibam", sample_out_bam,
               "-trackline", "-bga", "-5", "-strand", "-"]
        run_redirect_command(cmd, sample_out_prefix + "_minus.bedgraph")
        print("bedgraphs made on plus and minus strands")

    # Make length-sensitive alignments in H5 format.

    sample_out_h5 = sample_out_prefix + ".h5"
    second_id = config["SecondID"]
    if second_id is None:
        second_id = "NULL"
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "bam_to_h5.R"),
           "--Ncores=" + str(nprocesses),
           "--MinReadLen=" + str(config["MinReadLen"]),
           "--MaxReadLen=" + str(config["MaxReadLen"]),
           "--Buffer=" + str(config["Buffer"]),
           "--PrimaryID=" + config["PrimaryID"],
           "--SecondID=" + second_id,
           "--dataset=" + config["dataset"],
           "--bamFile=" + sample_out_bam,
           "--hdFile=" + sample_out_h5,
           "--orf_gff_file=" + config["orf_gff_file"],
           "--ribovizGFF=" + str(config["ribovizGFF"]),
           "--StopInCDS=" + str(config["StopInCDS"])]
    run_command(cmd)

    # Create summary statistics and analyses plots.

    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "generate_stats_figs.R"),
           "--Ncores=" + str(nprocesses),
           "--MinReadLen=" + str(config["MinReadLen"]),
           "--MaxReadLen=" + str(config["MaxReadLen"]),
           "--Buffer=" + str(config["Buffer"]),
           "--PrimaryID=" + config["PrimaryID"],
           "--dataset=" + config["dataset"],
           "--hdFile=" + sample_out_h5,
           "--out_prefix=" + sample_out_prefix,
           "--orf_fasta=" + config["orf_fasta"],
           "--rpf=" + str(config["rpf"]),
           "--orf_gff_file=" + config["orf_gff_file"],
           "--dir_out=" + out_dir,
           "--dir_data=" + data_dir,
           "--features_file=" + config["features_file"],
           "--do_pos_sp_nt_freq=" + str(config["do_pos_sp_nt_freq"])]
    run_command(cmd)
    print(("Finished processing sample " + fastq))


def collate_tpms(config_yaml, r_scripts):
    """
    Collate TPMs across sample results.

    :param config_yaml: YAML configuration file path
    :type config_yaml: str or unicode
    :param r_scripts:  Directory with RiboViz R scripts
    :type r_scripts: str or unicode
    :raise FileNotFoundError: if Rscript cannot be found
    :raise AssertionError: if collate_tpms.R returns non-zero exit
    code
    :raise FileNotFoundError: if a third-party tool cannot be found
    :raise AssertionError: if invocation of a third-party tool returns
    non-zero exit code
    """
    print("Collating TPMs across all processed samples")
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "collate_tpms.R"),
           "--yaml=" + config_yaml]
    run_command(cmd)


def prep_riboviz(py_scripts, r_scripts, data_dir, config_yaml):
    """
    Run the RiboViz workflow.

    Exit codes are as follows:

    * EXIT_OK (0): Processing successfully completed.
    * EXIT_CONFIG_ERROR (1): Errors occurred loading configuration.
    * EXIT_INDEX_ERROR (2): Error occurred during indexing.
    * EXIT_NO_SAMPLES_ERROR (3): No samples were provided.
    * EXIT_SAMPLES_ERROR (4): No sample was processed successfully.
    * EXIT_COLLATION_ERROR (5): Error occurred during TPMs collation.

    :param py_scripts: Directory with RiboViz Python scripts
    :type py_scripts: str or unicode
    :param r_scripts:  Directory with RiboViz R scripts
    :type r_scripts: str or unicode
    :param data_dir: Directory with data
    :type data_dir: str or unicode
    :param config_yaml: YAML configuration file path
    :type config_yaml: str or unicode
    :return: exit code
    :rtype: int
    """
    print(("Running under Python " + sys.version))
    # Extract configuration.
    try:
        with open(config_yaml, 'r') as f:
            config = yaml.load(f)
    except Exception:
        traceback.print_exc(file=sys.stdout)
        return EXIT_CONFIG_ERROR

    # Build indices for alignment, if necessary/requested.
    try:
        if config["build_indices"]:
            index_dir = config["dir_index"]
            if not os.path.exists(index_dir):
                os.makedirs(index_dir)
            r_rna_index = os.path.join(index_dir,
                                       config["rRNA_index"])
            build_indices(config["rRNA_fasta"], r_rna_index)
            print("rRNA index built")
            orf_index = os.path.join(index_dir,
                                     config["orf_index"])
            build_indices(config["orf_fasta"], orf_index)
            print("ORF index built")
    except KeyError:
        traceback.print_exc(file=sys.stdout)
        return EXIT_CONFIG_ERROR
    except Exception:
        traceback.print_exc(file=sys.stdout)
        return EXIT_INDEX_ERROR

    # Loop over sample fastq.gz files.
    print("Processing samples")
    try:
        in_dir = config["dir_in"]
        tmp_dir = config["dir_tmp"]
        out_dir = config["dir_out"]
    except KeyError:
        traceback.print_exc(file=sys.stdout)
        return EXIT_CONFIG_ERROR
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if ("fq_files" not in config) or \
       (config["fq_files"] is None) or \
       (len(config["fq_files"]) == 0):
        print("No samples are defined")
        return EXIT_NO_SAMPLES_ERROR
    samples = config["fq_files"]
    num_samples = len(config["fq_files"])
    failures = 0
    for sample in list(samples.keys()):
        try:
            fastq = os.path.join(in_dir, samples[sample])
            process_sample(sample,
                           fastq,
                           r_rna_index,
                           orf_index,
                           config,
                           py_scripts,
                           r_scripts,
                           data_dir,
                           tmp_dir,
                           out_dir)
        except Exception:
            traceback.print_exc(file=sys.stdout)
            failures += 1
    print(("Finished processing %d samples, %d failed"
           % (num_samples, failures)))
    if failures == num_samples:
        return EXIT_SAMPLES_ERROR

    # Collate TPMs across sample results.
    try:
        collate_tpms(config_yaml, r_scripts)
    except Exception:
        traceback.print_exc(file=sys.stdout)
        return EXIT_COLLATION_ERROR

    print("Completed")
    return EXIT_OK


if __name__ == "__main__":
    py_scripts = sys.argv[1]
    r_scripts = sys.argv[2]
    data_dir = sys.argv[3]
    config_yaml = sys.argv[4]
    exit_code = prep_riboviz(py_scripts, r_scripts, data_dir, config_yaml)
    sys.exit(exit_code)
