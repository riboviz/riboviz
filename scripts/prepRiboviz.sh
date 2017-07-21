#!/usr/bin/env bash

# prepRiboviz.sh: 
## Shell script for preparing ribosome profiling data for RiboViz or other analysis, that 
# - processes all fastq.gz files in an input directory  (-in), 
# - cuts out sequencing library adapters (CTGTAGGCACC or -adapters)
# - removes rRNA or other contaminating reads by hisat2 alignment to rRNA index file (-rRNA)
# - aligns remaining reads to ORFs or other hisat2 index file (-orf)
# - trims 5' mismatches from reads and removes reads with more than 2 mismatches
# - parallelizes over many processes (-p), except for cutadapt which isn't parallel
# - puts all intermediate files into a temporary directory (-tmp)
# - when finished, the script will put useful output files in another directory (-out)
# Note that the bamfiles ${fn_outbam} are directly usable in genome browsers, etc.

function usage
{
    echo "usage: prepRiboviz [-in dir_in] [-rRNA rRNA_index_file] [-orf orf_index_file]  [-out dir_out] [-a adapters] [-p processes] [-tmp tmp_dir]"
}

tmp_dir=tmp
adapters=CTGTAGGCACC
processes=4

#----- check for arguments -----
while [ "$1" != "" ]; do
    case $1 in
        -in )                   shift
                                dir_in=$1
                                ;;
        -rRNA )                 shift
                                rRNA_index_file=$1
                                ;;
        -orf )                  shift
                                orf_index_file=$1
                                ;;
        -out )                  shift
                                dir_out=$1
                                ;;
        -a )                    shift
                                adapters=$1
                                ;;
        -p )                    shift
                                processes=$1
                                ;;
        -tmp )                  shift
                                tmp_dir=$1
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

#----- print arguments -----
echo "prepRiboviz
rRNA_index_file = ${rRNA_index_file}
orf_index_file = ${orf_index_file}
dir_in = ${dir_in} 
dir_out = ${dir_out}
tmp_dir = ${tmp_dir}
adapters = ${adapters}"

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
mkdir ${tmp_dir}
mkdir ${dir_out}

#----- loop over fastq.gz files in directory -----
for fn in ${dir_in}/*.fastq.gz
    do
    ## report which file
    echo processing file ${fn}
    ## get filename stem
    fn_nodir=$(basename ${fn})
    fn_stem=${fn_nodir%%.fastq.gz}
    ## make filenames for:
    ##  trimmed reads
    fn_trim=${tmp_dir}/${fn_stem}_trim.fq
    ##  trimmed non-rRNA reads
    fn_nonrRNA=${tmp_dir}/${fn_stem}_nonrRNA.fq
    ##  rRNA-mapped reads
    fn_rRNA_mapped=${tmp_dir}/${fn_stem}_rRNA_map.sam
    ##  orf-mapped reads
    fn_orf_mapped=${tmp_dir}/${fn_stem}_orf_map.sam
    fn_orf_mapped_clean=${tmp_dir}/${fn_stem}_orf_map_clean.sam
    fn_nonaligned=${tmp_dir}/${fn_stem}_unaligned.sam
    ##  bam file
    fn_outbam=${dir_out}/${fn_stem}
    ##
    ## cut illumina adapters
    cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${fn_trim} ${fn}
    ## Map reads to rRNA
    echo hisat2 -p ${processes} -N 1 --un ${fn_nonrRNA} -x ${rRNA_index_file} \
    -S ${fn_rRNA_mapped} -U ${fn_trim}
    hisat2 -p ${processes} -N 1 --un ${fn_nonrRNA} -x ${rRNA_index_file} \
    -S ${fn_rRNA_mapped} -U ${fn_trim}
    # Map to orfs with (mostly) default settings, up to 2 alignments
    echo hisat2 -p ${processes} -k 2 --no-spliced-alignment --rna-strandness F --no-unal \
    --un ${fn_nonaligned} \
    -x ${orf_index_file} -S ${fn_orf_mapped} -U ${fn_nonrRNA}
    hisat2 -p ${processes} -k 2 --no-spliced-alignment --rna-strandness F --no-unal \
    --un ${fn_nonaligned} \
    -x ${orf_index_file} -S ${fn_orf_mapped} -U ${fn_nonrRNA}
    ## Trim 5' mismatched nt and remove reads with >1 mismatch
    python trim5pmismatch.py -mm 2 -in ${fn_orf_mapped} -out ${fn_orf_mapped_clean}
    ## convert sam (text) output to bam file (compressed binary)
    samtools view -b ${fn_orf_mapped_clean} | \
    ## sort bam file on genome and write
    samtools sort -@ ${processes} -O bam -o ${fn_outbam}.bam -
    ## index bamfile
    samtools index ${fn_outbam}.bam
    #####
    ## LOOK HERE!!!
    ## this is the point where we need to run reads_to_list, etc, to make length-sensitive alignments
    #####
    ## genome coverage as bedgraph
    ## calculate genome coverage for plus strand
    bedtools genomecov -ibam ${fn_outbam}.bam -bga -5 -strand + > ${fn_outbam}_plus.bedgraph
    ## calculate genome coverage for minus strand
    bedtools genomecov -ibam ${fn_outbam}.bam -bga -5 -strand - > ${fn_outbam}_minus.bedgraph
    ## ctrl-c interrupts entire loop, not just current iteration
    trap exit SIGINT
done

echo "prepRiboviz.sh completed"