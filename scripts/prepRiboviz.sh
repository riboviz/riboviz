#!/usr/bin/env bash

# prepRiboviz.sh: 
## Shell script for preparing ribosome profiling data for RiboViz or other analysis, that 
# - takes all configuration info from input yaml file.
# - processes all fastq.gz files named in config file (-in), 
# - cuts out sequencing library adapters (CTGTAGGCACC or -adapters)
# - removes rRNA or other contaminating reads by hisat2 alignment to rRNA index file (-rRNA)
# - aligns remaining reads to ORFs or other hisat2 index file (-orf)
# - trims 5' mismatches from reads and removes reads with more than 2 mismatches
# - parallelizes over many processes (-p), except for cutadapt which isn't parallel
# - puts all intermediate files into a temporary directory (-tmp)
# - when finished, the script will put useful output files in another directory (-out)
# Note that the bamfiles ${fn_outbam} are directly usable in genome browsers, etc.

# function usage
{
    echo "usage: prepRiboviz config.yaml"
}


# set directory name for script
dir_scripts=$(dirname $0)
# include parse_yaml function
. ${dir_scripts}/parse_yaml.sh
config_yaml=$1
config_yaml=vignette/vignette_config.yaml # just for debugging from command line

#----- print parameters from yaml -----
echo "prepRiboviz running with input parameters:"
parse_yaml ${config_yaml}
## here should check an input yaml file is supplied, with the right arguments.

#----- load parameters -----
eval $(parse_yaml ${config_yaml})

#----- build indices if necessary -----
if [ ${build_indices} == TRUE ]; then
    mkdir ${dir_index}
    hisat2-build ${rRNA_fasta} ${rRNA_index}
    hisat2-build ${orf_fasta} ${orf_index}
fi

#----- check there are sensible files in the right place -----
# if [ compgen -G "${dir_in}/*.fastq.gz" ];
if [ -z "ls ${dir_in}/*.fastq.gz" ];
then
    echo "directory ${dir_in} contains no fastq.qz files"
    exit 1
else 
    echo "fastq.gz files present"
fi


#----- create temp and output directory -----
mkdir ${dir_tmp}
mkdir ${dir_out}

#----- loop over fastq.gz files -----
fqfs=$(compgen -A variable | grep "fq_files")
for fq in ${fqfs}
    do
    echo processing fastq ${fq}
    ## get the filename from the variable
    eval fn_nodir=\$$fq
    fn=${dir_in}/${fn_nodir}
    ## report which file
    echo processing file ${fn_nodir}
    ## get filename stem
    # fn_nodir=$(basename ${fn})
    fn_stem=${fn_nodir%%.fastq.gz}
    ## make filenames for:
    ##  trimmed reads
    fn_trim=${dir_tmp}/${fn_stem}_trim.fq
    ##  trimmed non-rRNA reads
    fn_nonrRNA=${dir_tmp}/${fn_stem}_nonrRNA.fq
    ##  rRNA-mapped reads
    fn_rRNA_mapped=${dir_tmp}/${fn_stem}_rRNA_map.sam
    ##  orf-mapped reads
    fn_orf_mapped=${dir_tmp}/${fn_stem}_orf_map.sam
    fn_orf_mapped_clean=${dir_tmp}/${fn_stem}_orf_map_clean.sam
    fn_nonaligned=${dir_tmp}/${fn_stem}_unaligned.sam
    ##  bam + h5 files
    fn_out=${dir_out}/${fn_stem}
    ##
    ## cut illumina adapters
    cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${fn_trim} ${fn}
    ## Map reads to rRNA
    echo hisat2 -p ${nprocesses} -N 1 --un ${fn_nonrRNA} -x ${rRNA_index} \
    -S ${fn_rRNA_mapped} -U ${fn_trim}
    hisat2 -p ${nprocesses} -N 1 --un ${fn_nonrRNA} -x ${rRNA_index} \
    -S ${fn_rRNA_mapped} -U ${fn_trim}
    # Map to orfs with (mostly) default settings, up to 2 alignments
    echo hisat2 -p ${nprocesses} -k 2 --no-spliced-alignment --rna-strandness F --no-unal \
    --un ${fn_nonaligned} \
    -x ${orf_index} -S ${fn_orf_mapped} -U ${fn_nonrRNA}
    hisat2 -p ${nprocesses} -k 2 --no-spliced-alignment --rna-strandness F --no-unal \
    --un ${fn_nonaligned} \
    -x ${orf_index} -S ${fn_orf_mapped} -U ${fn_nonrRNA}
    ## Trim 5' mismatched nt and remove reads with >1 mismatch
    python ${dir_scripts}/trim_5p_mismatch.py -mm 2 -in ${fn_orf_mapped} -out ${fn_orf_mapped_clean}
    ## convert sam (text) output to bam file (compressed binary)
    samtools view -b ${fn_orf_mapped_clean} | \
    ## sort bam file on genome and write
    samtools sort -@ ${nprocesses} -O bam -o ${fn_out}.bam -
    ## index bamfile
    samtools index ${fn_out}.bam
    #####
    ## LOOK HERE!!!
    ## this is the point where we need to run reads_to_list, etc, to make length-sensitive alignments
    #     echo R CMD BATCH "--args Ncores=${nprocesses} SecondID=${SecondID} bamFile=${fn_outbam}.bam hdfile=${fn_outbam}.h5 gffFile=${orf_gff_file}" ${dir_scripts}/bam_to_h5.R
    #     R CMD BATCH "--args Ncores=${nprocesses} SecondID=${SecondID} bamFile=${fn_outbam}.bam hdfile=${fn_outbam}.h5 gffFile=${orf_gff_file}" ${dir_scripts}/bam_to_h5.R
    #  
    echo Rscript --vanilla ${dir_scripts}/bam_to_h5.R --Ncores=${nprocesses} --SecondID=${SecondID} --bamFile=${fn_outbam}.bam --hdFile=${fn_outbam}.h5 --orf_gff_file=${orf_gff_file}
    Rscript --vanilla ${dir_scripts}/bam_to_h5.R --Ncores=${nprocesses} --SecondID=${SecondID} --bamFile=${fn_outbam}.bam --hdFile=${fn_outbam}.h5 --orf_gff_file=${orf_gff_file}
    #####
    ## transcriptome coverage as bedgraph
    ## calculate transcriptome coverage for plus strand
    bedtools genomecov -ibam ${fn_out}.bam -bga -5 -strand + > ${fn_out}_plus.bedgraph
    ## calculate transcriptome coverage for minus strand
    bedtools genomecov -ibam ${fn_out}.bam -bga -5 -strand - > ${fn_out}_minus.bedgraph
    ## ctrl-c interrupts entire loop, not just current iteration
    trap exit SIGINT
done

echo "prepRiboviz.sh completed"