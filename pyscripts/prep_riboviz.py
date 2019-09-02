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
                   tmp_dir,
                   out_dir):
    """
    Process a single FASTQ sample file.

    :param sample: Sample name
    :type sample: str or unicode
    :param FASTQ: Sample FASTQ file
    :type FASTQ: str or unicode
    :param r_rna_index: Prefix of rRNA HT2 index files
    :type r_rna_index: str or unicode
    :param orf_index: Prefix of ORF HT2 index files
    :type orf_index: str or unicode
    :param config: RiboViz configuration
    :type config: dict
    :param tmp_dir Temporary files directory
    :type tmp_dir: str or unicode
    :param out_dir Output files directory
    :type out_dir: str or unicode
    """
    print(("Processing sample " + sample))
    # TODO raise and catch in caller.
    if not os.path.exists(fastq):
        print(("File " + fastq + " not found"))
        return
    print(("Processing file " + fastq))

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
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Map reads to rRNA.

    # Trimmed non-rRNA reads.
    non_r_rna_trim_fq = os.path.join(tmp_dir, sample + "_nonrRNA.fq")
    # rRNA-mapped reads.
    r_rna_map_sam = os.path.join(tmp_dir, sample + "_rRNA_map.sam")
    cmd = ["hisat2", "-p", str(config["nprocesses"]), "-N", "1",
           "--un", non_r_rna_trim_fq, "-x", r_rna_index,
           "-S", r_rna_map_sam, "-U", trim_fq]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Map to ORFs with (mostly) default settings, up to 2 alignments.

    # ORF-mapped reads.
    orf_map_sam = os.path.join(tmp_dir, sample + "_orf_map.sam")
    # Unaligned reads.
    unaligned_fq = os.path.join(tmp_dir, sample + "_unaligned.fq")
    cmd = ["hisat2", "-p", str(config["nprocesses"]), "-k", "2",
           "--no-spliced-alignment", "--rna-strandness",
           "F", "--no-unal", "--un", unaligned_fq,
           "-x", orf_index, "-S", orf_map_sam,
           "-U", non_r_rna_trim_fq]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Trim 5' mismatched nt and remove reads with >1 mismatch.

    # ORF-mapped reads.
    orf_map_sam_clean = os.path.join(tmp_dir,
                                     sample +
                                     "_orf_map_clean.sam")
    cmd = ["python", os.path.join(py_scripts, "trim_5p_mismatch.py"),
           "-mm", "2", "-in", orf_map_sam,
           "-out", orf_map_sam_clean]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Convert SAM (text) output to BAM (compressed binary).

    sample_out_prefix = os.path.join(out_dir, sample)
    sample_out_bam = sample_out_prefix + ".bam"
    # Use pattern suggested by following to implement pipe of
    # "samtools view" output into "samtools sort".
    # https://stackoverflow.com/questions/13332268/python-subprocess-command-with-pipe
    cmd_view = ["samtools", "view", "-b", orf_map_sam_clean]
    print(("Running: " + list_to_str(cmd_view)))
    # Sort BAM file on genome and write.
    cmd_sort = ["samtools", "sort", "-@", str(config["nprocesses"]),
                "-O", "bam", "-o", sample_out_bam, "-"]
    print(("Running: " + list_to_str(cmd_sort)))
    process_view = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
    output_sort = subprocess.check_output(cmd_sort, stdin=process_view.stdout)
    process_view.wait()

    # Index BAM file.
    cmd = ["samtools", "index", sample_out_bam]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    if config["make_bedgraph"]:
        # Record transcriptome coverage as a bedgraph.
        # Calculate transcriptome coverage for plus strand.
        # Use pattern suggested by following to implement stdout
        # redirection from "bedtools" to ".bedgraph" files.
        # https://www.saltycrane.com/blog/2008/09/how-get-stdout-and-stderr-using-python-subprocess-module/
        cmd = ["bedtools", "genomecov", "-ibam", sample_out_bam,
               "-trackline", "-bga", "-5", "-strand", "+"]
        print(("Running: " + list_to_str(cmd)))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, _ = p.communicate()
        with open(sample_out_prefix + "_plus.bedgraph", "wb") as f:
            f.write(stdout)
        # Calculate transcriptome coverage for minus strand.
        cmd = ["bedtools", "genomecov", "-ibam", sample_out_bam,
               "-trackline", "-bga", "-5", "-strand", "-"]
        print(("Running: " + list_to_str(cmd)))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, _ = p.communicate()
        with open(sample_out_prefix + "_minus.bedgraph", "wb") as f:
            f.write(stdout)
        print("bedgraphs made on plus and minus strands")

    # Make length-sensitive alignments in H5 format.

    sample_out_h5 = sample_out_prefix + ".h5"
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
           "--bamFile=" + sample_out_bam,
           "--hdFile=" + sample_out_h5,
           "--orf_gff_file=" + config["orf_gff_file"],
           "--ribovizGFF=" + str(config["ribovizGFF"]),
           "--StopInCDS=" + str(config["StopInCDS"])]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)

    # Create summary statistics and analyses plots.

    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "generate_stats_figs.R"),
           "--Ncores=" + str(config["nprocesses"]),
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
           "--orf_gff_file=" + config["orf_gff_file"],
           "--features_file=" + config["features_file"],
           "--do_pos_sp_nt_freq=" + str(config["do_pos_sp_nt_freq"])]
    print(("Running: " + list_to_str(cmd)))
    subprocess.call(cmd)
    print(("Finished processing sample " + fastq))

    
def collate_tpms(config_yaml):
    """
    Collate TPMs across sample results.

    :param config_yaml: YAML configuration file path
    :type config_yaml: str or unicode
    :raise FileNotFoundError: if Rscript cannot be found
    :raise AssertionError: if collate_tpms.R returns non-zero exit
    code
    """
    print("Collating TPMs across all processed amples")
    cmd = ["Rscript", "--vanilla",
           os.path.join(r_scripts, "collate_tpms.R"),
           "--yaml=" + config_yaml]
    exit_code = subprocess.call(cmd)
    assert exit_code == 0, "%s failed with exit code %d" % (cmd, exit_code)


def prep_riboviz(py_scripts, r_scripts, data_dir, config_yaml):
    """
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
    with open(config_yaml, 'r') as f:
        config = yaml.load(f)

    # Build indices for alignment, if necessary/requested.
    if config["build_indices"]:
        if not os.path.exists(config["dir_index"]):
            os.makedirs(config["dir_index"])
        r_rna_index = os.path.join(config["dir_index"], config["rRNA_index"])
        build_indices(config["rRNA_fasta"], r_rna_index)
        print("rRNA index built")
        orf_index = os.path.join(config["dir_index"], config["orf_index"])
        build_indices(config["orf_fasta"], orf_index)
        print("orf index built")

    if len(glob.glob(os.path.join(config["dir_in"], "*.fastq.gz"))) == 0:
        print(("Directory " +
               config["dir_in"] +
               " contains no fastq.qz files"))
        exit(1)

    if not os.path.exists(config["dir_tmp"]):
        os.makedirs(config["dir_tmp"])
    if not os.path.exists(config["dir_out"]):
        os.makedirs(config["dir_out"])

    # Loop over sample fastq.gz files.
    print("Processing samples")
    for sample in list(config["fq_files"].keys()):
        fastq = os.path.join(config["dir_in"], config["fq_files"][sample])
        process_sample(sample,
                       fastq,
                       r_rna_index,
                       orf_index,
                       config,
                       config["dir_tmp"],
                       config["dir_out"])
    print("Finished processing samples")

    # Collate TPMs across sample results.
    collate_tpms(config_yaml)

    print("Completed")
    return 0


if __name__ == "__main__":
    py_scripts = sys.argv[1]
    r_scripts = sys.argv[2]
    data_dir = sys.argv[3]
    config_yaml = sys.argv[4]
    prep_riboviz(py_scripts, r_scripts, data_dir, config_yaml)
