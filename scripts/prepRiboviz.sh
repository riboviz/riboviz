#!/usr/bin/env bash

# prepRiboviz.sh: 
## Shell script for preparing ribosome profiling data for RiboViz or other analysis, that 
# - takes all configuration info from input yaml file.
# - builds hisat2 indices if requested (build_indices is TRUE), in index directory (dir_index).
# - processes all fastq.gz files named in config file (dir_in), 
# - cuts out sequencing library adapters (CTGTAGGCACC or adapters)
# - removes rRNA or other contaminating reads by hisat2 alignment to rRNA index file (rRNA_index)
# - aligns remaining reads to ORFs or other hisat2 index file (orf_index)
# - trims 5' mismatches from reads and removes reads with more than 2 mismatches
# - parallelizes over many processes (nprocesses), except for cutadapt which isn't parallel
# - makes length-sensitive alignments in compressed h5 format by running reads\_to\_list.R
# - generates summary statistics, and analyses and QC plots for both RPF and mRNA datasets, by running generate\_stats\_figs.R
# - puts all intermediate files into a temporary directory (dir_tmp)
# - when finished, the script will put useful output files in another directory (dir_out)
# Note that the bamfiles ${fn_outbam} are directly usable in genome browsers, etc.
# Also optionally exports bedgraph files for plus and minus strands, if (make_bedgraph) is TRUE

# function usage
{
    echo "usage: prepRiboviz config.yaml"
}


# set directory name for script
dir_scripts=$(dirname $0)
# include parse_yaml function
. ${dir_scripts}/parse_yaml.sh
config_yaml=$1

## for debugging from command line staring at Riboviz directory
# dir_scripts=scripts
# config_yaml=vignette/vignette_config.yaml 

#----- print parameters from yaml -----
echo "prepRiboviz running with input parameters:"
parse_yaml ${config_yaml}
## here should check an input yaml file is supplied, with the right arguments.

#----- load parameters -----
eval $(parse_yaml ${config_yaml})

#----- build indices if necessary -----
if [ ${build_indices} == TRUE ]; then
    echo 
    echo "Building indices for alignment"
    mkdir -p ${dir_index} # this is not really used; needs fixing
    echo hisat2-build ${rRNA_fasta} ${rRNA_index}
    hisat2-build ${rRNA_fasta} ${rRNA_index}
    echo "rRNA index built"
    echo
    echo hisat2-build ${orf_fasta} ${orf_index}
    hisat2-build ${orf_fasta} ${orf_index}
    echo "orf index built"
    echo
fi

#----- check there are sensible files in the right place -----
# if [ compgen -G "${dir_in}/*.fastq.gz" ];
if [ -z "ls ${dir_in}/*.fastq.gz" ];
then
    echo "directory ${dir_in} contains no fastq.qz files"
    exit 1
else 
    echo "fastq.gz files present in ${dir_in}"
    echo
fi


#----- create temp and output directory -----
mkdir -p ${dir_tmp}
mkdir -p ${dir_out}

#----- loop over fastq.gz files -----
fqfs=$(compgen -A variable | grep "fq_files")
for fq in ${fqfs}
    do  
    echo processing fastq sample ${fq}
    ## get the filename from the variable
    eval fn_nodir=\$$fq
    fn=${dir_in}/${fn_nodir}
    ## check file present
    if [ ! -f ${fn} ]; then 
        ## break loop if absent
        echo File ${fn} not found
        echo 
        continue
    else  
        ## report which file
        echo processing file ${fn_nodir} 
        ## get filename stem
        # fn_stem=${fn_nodir%%.fastq.gz} # use fastq file names as file prefix
        eval fn_stem=${fq##fq_files_}	# use user-defined dataset name as file prefix
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
        cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${fn_trim} ${fn} -j ${nprocesses}
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
        echo python ${dir_scripts}/trim_5p_mismatch.py -mm 2 -in ${fn_orf_mapped} -out ${fn_orf_mapped_clean}
        python ${dir_scripts}/trim_5p_mismatch.py -mm 2 -in ${fn_orf_mapped} -out ${fn_orf_mapped_clean}
        ## convert sam (text) output to bam file (compressed binary)
        echo samtools view -b ${fn_orf_mapped_clean} 
        echo samtools sort -@ ${nprocesses} -O bam -o ${fn_out}.bam -
        samtools view -b ${fn_orf_mapped_clean} | \
        ## sort bam file on genome and write
        samtools sort -@ ${nprocesses} -O bam -o ${fn_out}.bam -
        ## index bamfile
        samtools index ${fn_out}.bam
        ##
        if [ ${make_bedgraph} == TRUE ]; then
            ## transcriptome coverage as bedgraph
            ## calculate transcriptome coverage for plus strand
            bedtools genomecov -ibam ${fn_out}.bam -bga -5 -strand + > ${fn_out}_plus.bedgraph
            ## calculate transcriptome coverage for minus strand
            bedtools genomecov -ibam ${fn_out}.bam -bga -5 -strand - > ${fn_out}_minus.bedgraph
            echo bedgraphs made on plus and minus strands
        fi
        ## run reads_to_list to make length-sensitive alignments in h5 format
        echo Rscript --vanilla ${dir_scripts}/bam_to_h5.R \
        --Ncores=${nprocesses} --MinReadLen=${MinReadLen} --MaxReadLen=${MaxReadLen} --Buffer=${Buffer} \
        --PrimaryID=${PrimaryID} --SecondID=${SecondID} --dataset=${dataset} --bamFile=${fn_out}.bam --hdFile=${fn_out}.h5 \
        --orf_gff_file=${orf_gff_file} --ribovizGFF=${ribovizGFF} --StopInCDS=${StopInCDS}
        #
        Rscript --vanilla ${dir_scripts}/bam_to_h5.R \
        --Ncores=${nprocesses} --MinReadLen=${MinReadLen} --MaxReadLen=${MaxReadLen} --Buffer=${Buffer} \
        --PrimaryID=${PrimaryID} --SecondID=${SecondID} --dataset=${dataset} --bamFile=${fn_out}.bam --hdFile=${fn_out}.h5 \
        --orf_gff_file=${orf_gff_file} --ribovizGFF=${ribovizGFF} --StopInCDS=${StopInCDS}
        ##
        ## generate summary statistics and analyses plots
        echo Rscript --vanilla ${dir_scripts}/generate_stats_figs.R \
        --Ncores=${nprocesses} --MinReadLen=${MinReadLen} --MaxReadLen=${MaxReadLen} --Buffer=${Buffer} \
        --PrimaryID=${PrimaryID} --dataset=${dataset} --hdFile=${fn_out}.h5 --out_prefix=${fn_out} --orf_fasta=${orf_fasta} \
        --orf_gff_file=${orf_gff_file} --rpf=${rpf} --dir_out=${dir_out} --dir_scripts=${dir_scripts} \
        --features_file=${features_file} --do_pos_sp_nt_freq=${do_pos_sp_nt_freq}
        #
        Rscript --vanilla ${dir_scripts}/generate_stats_figs.R \
        --Ncores=${nprocesses} --MinReadLen=${MinReadLen} --MaxReadLen=${MaxReadLen} --Buffer=${Buffer} \
        --PrimaryID=${PrimaryID} --dataset=${dataset} --hdFile=${fn_out}.h5 --out_prefix=${fn_out} --orf_fasta=${orf_fasta} \
        --orf_gff_file=${orf_gff_file} --rpf=${rpf} --dir_out=${dir_out} --dir_scripts=${dir_scripts} \
        --features_file=${features_file}  --do_pos_sp_nt_freq=${do_pos_sp_nt_freq}
        ##
        echo finished processing sample ${fq}
        echo
    fi
    ## ctrl-c interrupts entire loop, not just current iteration
    trap exit SIGINT
done

echo collating TPMs across samples
echo Rscript --vanilla ${dir_scripts}/collate_tpms.R --yaml=${config_yaml}
Rscript --vanilla ${dir_scripts}/collate_tpms.R --yaml=${config_yaml}

echo "prepRiboviz.sh completed"
echo
