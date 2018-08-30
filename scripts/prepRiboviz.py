#!/usr/bin/env python
"""
RiboViz workflow.

Usage:

    python scripts/prepRiboViz.py <SCRIPTS_DIRECTORY> <YAML_CONFIG_FILE>

Example:

    python scripts/prepRiboviz.py scripts/ vignette/vignette_config.yaml

Prepare ribosome profiling data for RiboViz or other analysis:

* - Read configuration info from YAML file.
* - Build hisat2 indices if requested (config["build_indices"] == True), in index directory (config["dir_index"]).
* - Process all fastq.gz files (config["dir_in"]).
* - Cut out sequencing library adapters ("CTGTAGGCACC" or config["adapters").
* - Remove rRNA or other contaminating reads by hisat2 alignment to rRNA index file (config["rRNA_index"]).
* - Align remaining reads to ORFs or other hisat2 index file (config["orf_index"])
* - Trim 5' mismatches from reads and remove reads with more than 2 mismatches
* - Parallelize over many processes (config["nprocesses"]), except for cutadapt which isn't parallel.
* - Make length-sensitive alignments in compressed h5 format by running "reads_to_list.R".
* - Generate summary statistics, and analyses and QC plots for both RPF and mRNA datasets, by running "generate_stats_figs.R".
* - Put all intermediate files into a temporary directory (config["dir_tmp"]).
* - When finished, put useful output files into output directory (config["dir_out"]).
* Optionally export bedgraph files for plus and minus strands, if requested (config["make_bedgraph"] == True).

bamfiles (config["dir_out"]/*.bam) are directly usable in genome browsers, etc.
"""

# RPF sequences => .fastq => trim & QC => .fastq => rRNA removed => .fastq -> transcriptome alignment => .bam => count 5' ends => .bed => bam_to_h5.R => .h5 => summary stats => R Shiny, PDF, etc

import glob
import os
import os.path
import subprocess
import sys
import yaml


def list_to_str(lst):
    return ' '.join(map(str, lst))


dir_scripts = sys.argv[1]
config_yaml = sys.argv[2] # vignette/vignette_config.yaml
with open(config_yaml, 'r') as f:
    config = yaml.load(f)

# Build indices if necessary.
if config["build_indices"]:
    # Build indices for alignment
    if not os.path.exists(config["dir_index"]): # "vignette/index"
        os.makedirs(config["dir_index"])
    # hisat2-build vignette/input/yeast_rRNA_R64-1-1.fa vignette/index/yeast_rRNA
    cmd = ["hisat2-build", config["rRNA_fasta"], config["rRNA_index"]]
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)
    print("rRNA index built")

    # hisat2-build vignette/input/yeast_YAL_CDS_w_250utrs.fa vignette/index/YAL_CDS_w_250
    cmd = ["hisat2-build", config["orf_fasta"], config["orf_index"]]
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)
    print("orf index built")

if len(glob.glob(os.path.join(config["dir_in"], "*.fastq.gz")))== 0:
    # config["dir_in"] "vignette/input"
    print("Directory " + config["dir_in"] + " contains no fastq.qz files")
    exit(1)

if not os.path.exists(config["dir_tmp"]): # "vignette/tmp"
    os.makedirs(config["dir_tmp"])
if not os.path.exists(config["dir_out"]): # "vignette/output"
    os.makedirs(config["dir_out"])

# Loop over fastq.gz files.
for fq_file in config["fq_files"].keys(): # [WTnone, WT3AT, NotHere]
    print("Processing fastq sample " + fq_file)
    # Get file name.
    fn_nodir = config["fq_files"][fq_file] # SRR1042855_s1mi.fastq.gz, SRR1042864_s1mi.fastq.gz, example_missing_file.fastq.gz
    fn = os.path.join(config["dir_in"], fn_nodir) # vignette/input/SRR1042855_s1mi.fastq.gz,...
    if not os.path.exists(fn):
        print("File " + fn + " not found")
        continue
    print("Processing file " + fn)

    # Use user-defined dataset name as file prefix.
    fn_stem = fq_file # WTnone,...
    # Create tmp and out file names...
    # Trimmed reads.
    fn_trim = os.path.join(config["dir_tmp"], fn_stem + "_trim.fq") # vignette/tmp/WTnone_trim.fq
    # Trimmed non-rRNA reads.
    fn_nonrRNA = os.path.join(config["dir_tmp"], fn_stem + "_nonrRNA.fq") # vignette/tmp/WTnone_nonrRNA.fq
    # rRNA-mapped reads.
    fn_rRNA_mapped = os.path.join(config["dir_tmp"], fn_stem + "_rRNA_map.sam") # vignette/tmp/WTnone_rRNA_map.sam
    # orf-mapped reads.
    fn_orf_mapped = os.path.join(config["dir_tmp"], fn_stem + "_orf_map.sam") # vignette/tmp/WTnone_orf_map.sam
    fn_orf_mapped_clean = os.path.join(config["dir_tmp"], fn_stem + "_orf_map_clean.sam") # vignette/tmp/WTnone_orf_map_clean.sam
    fn_nonaligned = os.path.join(config["dir_tmp"], fn_stem + "_unaligned.sam") # vignette/tmp/WTnone_unaligned.sam
    # bam and h5 files.
    fn_out = os.path.join(config["dir_out"], fn_stem) # vignette/output/WTnone

    # Cut illumina adapters.
    # cutadapt --trim-n -O 1 -m 5 -a CTGTAGGCACC -o vignette/tmp/WTnone_trim.fq vignette/input/SRR1042855_s1mi.fastq.gz -j 1
    cmd = ["cutadapt", "--trim-n", "-O", "1", "-m", "5", "-a", config["adapters"], "-o", fn_trim, fn, "-j", str(config["nprocesses"])]
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)

    # Map reads to rRNA.
    # hisat2 -p 1 -N 1 --un vignette/tmp/WTnone_nonrRNA.fq -x vignette/index/yeast_rRNA -S vignette/tmp/WTnone_rRNA_map.sam -U vignette/tmp/WTnone_trim.fq
    cmd = ["hisat2", "-p", str(config["nprocesses"]), "-N", "1", "--un", fn_nonrRNA, "-x", config["rRNA_index"], "-S", fn_rRNA_mapped, "-U", fn_trim]
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)

    # Map to orfs with (mostly) default settings, up to 2 alignments.
    # hisat2 -p 1 -k 2 --no-spliced-alignment --rna-strandness F --no-unal --un vignette/tmp/WTnone_unaligned.sam -x vignette/index/YAL_CDS_w_250 -S vignette/tmp/WTnone_orf_map.sam -U vignette/tmp/WTnone_nonrRNA.fq
    cmd = ["hisat2", "-p", str(config["nprocesses"]), "-k", "2", "--no-spliced-alignment", "--rna-strandness", "F", "--no-unal", "--un", fn_nonaligned, "-x", config["orf_index"], "-S", fn_orf_mapped, "-U", fn_nonrRNA]
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)

    # Trim 5' mismatched nt and remove reads with >1 mismatch.
    # python scripts/trim_5p_mismatch.py -mm 2 -in vignette/tmp/WTnone_orf_map.sam -out vignette/tmp/WTnone_orf_map_clean.sam
    cmd = ["python", os.path.join(dir_scripts, "trim_5p_mismatch.py"), "-mm", "2", "-in", fn_orf_mapped, "-out", fn_orf_mapped_clean]
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)

    # Convert sam (text) output to bam file (compressed binary).
    # samtools view -b vignette/tmp/WTnone_orf_map_clean.sam | samtools sort -@ 1 -O bam -o vignette/output/WTnone.bam -
    # Use pattern suggested by following to implement pipe:
    # https://stackoverflow.com/questions/13332268/python-subprocess-command-with-pipe
    cmd_view = ["samtools", "view", "-b", fn_orf_mapped_clean]
    print("Running: " + list_to_str(cmd_view))
    cmd_sort = ["samtools", "sort", "-@", str(config["nprocesses"]), "-O", "bam", "-o", fn_out + ".bam", "-"]
    print("Running: " + list_to_str(cmd_sort))
    process_view = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
    output_sort = subprocess.check_output(cmd_sort, stdin=process_view.stdout)
    process_view.wait()

    # Sort bam file on genome and write.
    # Index bamfile.
    # samtools index vignette/output/WTnone.bam
    cmd = ["samtools", "index", fn_out + ".bam"]
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)

    if config["make_bedgraph"]:
        # Transcriptome coverage as bedgraph.
        # Calculate transcriptome coverage for plus strand.
        # cmd = ["bedtools", "genomecov", "-ibam", fn_out + ".bam", "-bga", "-5", "-strand", "+", ">", fn_out + "_plus.bedgraph"]
        # Use pattern suggested by following to implement stdout redirection:
        # https://www.saltycrane.com/blog/2008/09/how-get-stdout-and-stderr-using-python-subprocess-module/
        cmd = ["bedtools", "genomecov", "-ibam", fn_out + ".bam", "-bga", "-5", "-strand", "+"]
        print("Running: " + list_to_str(cmd))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, _ = p.communicate()
        with open(fn_out + "_plus.bedgraph", "w") as f:
            f.write(stdout)
        # Calculate transcriptome coverage for minus strand.
        # cmd = ["bedtools", "genomecov", "-ibam", fn_out + ".bam", "-bga", "-5", "-strand", "-", ">", fn_out + "_minus.bedgraph"]
        cmd = ["bedtools", "genomecov", "-ibam", fn_out + ".bam", "-bga", "-5", "-strand", "-"]
        print("Running: " + list_to_str(cmd))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, _ = p.communicate()
        with open(fn_out + "_minus.bedgraph", "w") as f:
            f.write(stdout)
    print("bedgraphs made on plus and minus strands")

    # Run reads_to_list to make length-sensitive alignments in h5 format.
    # Rscript --vanilla scripts/bam_to_h5.R --Ncores=1 --MinReadLen=10 --MaxReadLen=50 --Buffer=250 --PrimaryID=Name --SecondID=NULL --dataset=vignette --bamFile=vignette/output/WTnone.bam --hdFile=vignette/output/WTnone.h5 --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 --ribovizGFF=TRUE --StopInCDS=FALSE
    second_id = config["SecondID"]
    if second_id is None:
        second_id = "NULL"
    cmd = ["Rscript", "--vanilla",
        os.path.join(dir_scripts, "bam_to_h5.R"),
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
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)
    # Generate summary statistics and analyses plots.
    # Rscript --vanilla scripts/generate_stats_figs.R --Ncores=1 --MinReadLen=10 --MaxReadLen=50 --Buffer=250 --PrimaryID=Name --dataset=vignette --hdFile=vignette/output/WTnone.h5 --out_prefix=vignette/output/WTnone --orf_fasta=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=TRUE --dir_out=vignette/output --dir_scripts=scripts
    cmd = ["Rscript", "--vanilla",
        os.path.join(dir_scripts, "generate_stats_figs.R"),
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
        "--dir_out=" + config["dir_out"],
        "--dir_scripts=" + dir_scripts]
    print("Running: " + list_to_str(cmd))
    subprocess.call(cmd)
    print("Finished processing sample " + fq_file)
