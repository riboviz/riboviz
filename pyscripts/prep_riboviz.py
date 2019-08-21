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
*   rRNA index file (config["rRNA_index"]).
* - Align remaining reads to ORFs or other hisat2 index file
*   (config["orf_index"]).
* - Trim 5' mismatches from reads and remove reads with more than 2
*   mismatches.
* - Parallelize over many processes (config["nprocesses"]), except for
*   cutadapt which isn't parallel.
* - Make length-sensitive alignments in compressed h5 format by
*   running "bam_to_h5.R".
* - Generate summary statistics, and analyses and QC plots for both
*   RPF and mRNA datasets, by running "generate_stats_figs.R".
* - Put all intermediate files into a temporary directory
*   (config["dir_tmp"]).
* - When finished, put useful output files into output directory
*   (config["dir_out"]).
*
* Optionally export bedgraph files for plus and minus strands, if
*   requested (config["make_bedgraph"] == True).

bamfiles (config["dir_out"]/*.bam) are directly usable in genome browsers, etc.
"""

import glob
import os
import os.path
import subprocess
import sys
import yaml


def list_to_str(lst):
    """
    Convert list to space-delimited string.

    :param lst: list
    :type lst: list
    :return: list as string
    :rtype: str or unicode
    """
    return ' '.join(map(str, lst))


py_scripts = sys.argv[1]
r_scripts = sys.argv[2]
data_dir = sys.argv[3]
config_yaml = sys.argv[4]
with open(config_yaml, 'r') as f:
    config = yaml.load(f)

# Build indices if necessary.
if config["build_indices"]:
    # Build indices for alignment
    if not os.path.exists(config["dir_index"]):
        os.makedirs(config["dir_index"])

    cmd = ["hisat2-build", config["rRNA_fasta"], config["rRNA_index"]]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)
    print("rRNA index built")

    cmd = ["hisat2-build", config["orf_fasta"], config["orf_index"]]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)
    print("orf index built")

if len(glob.glob(os.path.join(config["dir_in"], "*.fastq.gz"))) == 0:
    print(("Directory " + config["dir_in"] + " contains no fastq.qz files"))
    exit(1)

if not os.path.exists(config["dir_tmp"]):
    os.makedirs(config["dir_tmp"])
if not os.path.exists(config["dir_out"]):
    os.makedirs(config["dir_out"])

# Loop over fastq.gz files.
for fq_file in list(config["fq_files"].keys()):
    print(("Processing fastq sample " + fq_file))
    # Get file name.
    fn_nodir = config["fq_files"][fq_file]
    fn = os.path.join(config["dir_in"], fn_nodir)
    if not os.path.exists(fn):
        print(("File " + fn + " not found"))
        continue
    print(("Processing file " + fn))

    # Use user-defined dataset name as file prefix.
    fn_stem = fq_file
    # Create tmp and out file names...
    # Trimmed reads.
    fn_trim = os.path.join(config["dir_tmp"],
                           fn_stem + "_trim.fq")
    # Trimmed non-rRNA reads.
    fn_nonrRNA = os.path.join(config["dir_tmp"],
                              fn_stem + "_nonrRNA.fq")
    # rRNA-mapped reads.
    fn_rRNA_mapped = os.path.join(config["dir_tmp"],
                                  fn_stem + "_rRNA_map.sam")
    # orf-mapped reads.
    fn_orf_mapped = os.path.join(config["dir_tmp"],
                                 fn_stem + "_orf_map.sam")
    fn_orf_mapped_clean = os.path.join(config["dir_tmp"],
                                       fn_stem + "_orf_map_clean.sam")
    fn_nonaligned = os.path.join(config["dir_tmp"],
                                 fn_stem + "_unaligned.fq")
    # bam and h5 files.
    fn_out = os.path.join(config["dir_out"], fn_stem)

    # Cut illumina adapters.
    cmd = ["cutadapt", "--trim-n", "-O", "1", "-m", "5",
           "-a", config["adapters"], "-o", fn_trim, fn,
           "-j", str(config["nprocesses"])]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Map reads to rRNA.
    cmd = ["hisat2", "-p", str(config["nprocesses"]), "-N", "1",
           "--un", fn_nonrRNA, "-x", config["rRNA_index"],
           "-S", fn_rRNA_mapped, "-U", fn_trim]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Map to orfs with (mostly) default settings, up to 2 alignments.
    cmd = ["hisat2", "-p", str(config["nprocesses"]), "-k", "2",
           "--no-spliced-alignment", "--rna-strandness",
           "F", "--no-unal", "--un", fn_nonaligned,
           "-x", config["orf_index"], "-S", fn_orf_mapped,
           "-U", fn_nonrRNA]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Trim 5' mismatched nt and remove reads with >1 mismatch.
    cmd = ["python", os.path.join(py_scripts, "trim_5p_mismatch.py"),
           "-mm", "2", "-in", fn_orf_mapped,
           "-out", fn_orf_mapped_clean]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Convert sam (text) output to bam file (compressed binary).
    # Use pattern suggested by following to implement pipe of
    # "samtools view" output into "samtools sort".
    # https://stackoverflow.com/questions/13332268/python-subprocess-command-with-pipe
    cmd_view = ["samtools", "view", "-b", fn_orf_mapped_clean]
    print(("Running: " + list_to_str(cmd_view)))
    # Sort bam file on genome and write.
    cmd_sort = ["samtools", "sort", "-@", str(config["nprocesses"]),
                "-O", "bam", "-o", fn_out + ".bam", "-"]
    print(("Running: " + list_to_str(cmd_sort)))
    process_view = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
    output_sort = subprocess.check_output(cmd_sort, stdin=process_view.stdout)
    process_view.wait()

    # Index bamfile.
    cmd = ["samtools", "index", fn_out + ".bam"]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    if config["make_bedgraph"]:
        # Transcriptome coverage as bedgraph.
        # Calculate transcriptome coverage for plus strand.
        # Use pattern suggested by following to implement stdout
        # redirection from "bedtools" to ".bedgraph" files.
        # https://www.saltycrane.com/blog/2008/09/how-get-stdout-and-stderr-using-python-subprocess-module/
        cmd = ["bedtools", "genomecov", "-ibam", fn_out + ".bam",
               "-trackline", "-bga", "-5", "-strand", "+"]
        print(("Running: " + list_to_str(cmd)))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, _ = p.communicate()
        with open(fn_out + "_plus.bedgraph", "wb") as f:
            f.write(stdout)
        # Calculate transcriptome coverage for minus strand.
        cmd = ["bedtools", "genomecov", "-ibam", fn_out + ".bam",
               "-trackline", "-bga", "-5", "-strand", "-"]
        print(("Running: " + list_to_str(cmd)))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, _ = p.communicate()
        with open(fn_out + "_minus.bedgraph", "wb") as f:
            f.write(stdout)
    print("bedgraphs made on plus and minus strands")

    # Run bam_to_h5 to make length-sensitive alignments in h5 format.
    second_id = config["SecondID"]
    if second_id is None:
        second_id = "NULL"
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "bam_to_h5.R"),
           "--Ncores=" + str(config["nprocesses"]),
           "--MinReadLen=" + str(config["MinReadLen"]),
           "--MaxReadLen=" + str(config["MaxReadLen"]),
           "--Buffer=" + str(config["Buffer"]),
           "--PrimaryID=" + config["PrimaryID"],
           "--SecondID=" + second_id,
           "--dataset=" + config["dataset"],
           "--bamFile=" + fn_out + ".bam",
           "--hdFile=" + fn_out + ".h5",
           "--orf_gff_file=" + config["orf_gff_file"],
           "--ribovizGFF=" + str(config["ribovizGFF"]),
           "--StopInCDS=" + str(config["StopInCDS"])]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)
    # Generate summary statistics and analyses plots.
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "generate_stats_figs.R"),
           "--Ncores=" + str(config["nprocesses"]),
           "--MinReadLen=" + str(config["MinReadLen"]),
           "--MaxReadLen=" + str(config["MaxReadLen"]),
           "--Buffer=" + str(config["Buffer"]),
           "--PrimaryID=" + config["PrimaryID"],
           "--dataset=" + config["dataset"],
           "--hdFile=" + fn_out + ".h5",
           "--out_prefix=" + fn_out,
           "--orf_fasta=" + config["orf_fasta"],
           "--rpf=" + str(config["rpf"]),
           "--orf_gff_file=" + config["orf_gff_file"],
           "--dir_out=" + config["dir_out"],
           "--dir_data=" + data_dir,
           "--orf_gff_file=" + config["orf_gff_file"],
           "--features_file=" + config["features_file"],
           "--do_pos_sp_nt_freq=" + str(config["do_pos_sp_nt_freq"])]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)
    print(("Finished processing sample " + fq_file))

print("collating TPMs across samples")
cmd = ["Rscript", "--vanilla",
	   os.path.join(r_scripts, "collate_tpms.R"),
	   "--yaml=" + config_yaml ]
subprocess.call(cmd)

print("finished running prepRiboviz.py")
