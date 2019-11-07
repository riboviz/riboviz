#! python

## SplitRNATagSeqFastQ.py
## Assigns RNATagSeq reads to samples based on initial 8-nt barcode (TagRead)
## Inputs:
##  - One .fastq.gz file with sequencing reads, 
##  OR - Two paired .fastq.gz files with read pairs in corresponding positions
##  - r1 file should have TagReads at beginning (5'/left)
##  - tab delimited .txt file with columns (at least) SampleID, TagRead
##  - in principle this should harvest the P7 read information from fastq header; but it doesn't
##
## Outputs:
##  - for each SampleID, two paired .fastq.gz files containing assigned reads, ending R1 and R2
## 
## Edward Wallace ewjwallace@gmail.com, 2017

import sys, os, csv, gzip, shutil, argparse
from itertools import islice
import pandas as pd
import re 

def trim_fastq_record(fqr,n=9):
    """trim initial n letters from fastq record"""
    return [ fqr[0], fqr[1][n:], fqr[2], fqr[3][n:] ]

def startswith_mismatch_regex(string,barcode,mm=0):
    """return true if fastq record starts with barcode, up to mm mismatches"""
    # same input-output as startswith_mismatch, but uses regex. slower.
    teststring = "\G(" + barcode + "){s<=" + str(mm) + "}"
    if re.findall(teststring, string) :
        return True
    else :
        return False

def startswith_mismatch(string,barcode,mm=0):
    """return true if string record starts with barcode, up to mm mismatches"""
    m = sum([string[i] != barcode[i] for i in range(len(barcode))])
    return m <= mm


if __name__=="__main__" :    
    # test1sing: python SplitRNATagSeqFastQ.py -r1 data/Sample_4reads_R1.fastq.gz -ss data/TagSeqBarcodedOligos2015.txt -o TestSingleSplit4reads
    # test2sing: python SplitRNATagSeqFastQ.py -r1 data/Sample_init10000_R1.fastq.gz -ss data/TagSeqBarcodedOligos2015.txt -o TestSingleSplit10000
    # test1pair: python SplitRNATagSeqFastQ.py -r1 data/Sample_4reads_R1.fastq.gz -r2 data/Sample_4reads_R2.fastq.gz -ss data/TagSeqBarcodedOligos2015.txt -o TestPairSplit4reads
    # test2pair: python SplitRNATagSeqFastQ.py -r1 data/Sample_init10000_R1.fastq.gz -r2 data/Sample_init10000_R2.fastq.gz -ss data/TagSeqBarcodedOligos2015.txt -o TestPairSplit10000
        
    # define input options
    parser = argparse.ArgumentParser(description="Demultiplex reads from fastq.gz by inline barcodes")
    parser.add_argument("-r1", "--read1",dest="r1_fn",nargs='?',help="read 1 filename in fastq.gz format")
    parser.add_argument("-r2", "--read2",dest="r2_fn",default=False,nargs='?',help="read 2 pair filename in fastq.gz format")
    parser.add_argument("-to", "--trimout",dest="trimout",default=True,nargs='?',help="trim initial TagRead from read 1")
    parser.add_argument("-ss", "--samplesheet", dest="samplesheet_fn",nargs='?',help="samplesheet filename, tab-delim text format with SampleID and TagRead columns")
    parser.add_argument("-o", "--outdir", dest="outdir",nargs='?',default="TestSplitRNATagSeq",help="output directory")
    parser.add_argument("-m", "--mismatches", dest="mismatches", default=1, type=int, help="number of mismatches permitted in barcode")
    options = parser.parse_args()
    
    print "Demultiplexing reads for file:\n" + options.r1_fn + \
        "\nusing sample sheet:\n" + options.samplesheet_fn
    
    if not os.path.isfile(options.samplesheet_fn):
        raise IOError("# Error: sample sheet file {} does not exist".format(options.samplesheet_fn))
    # read sample sheet in using pandas package
    samplesheet = pd.read_csv(options.samplesheet_fn,comment="#",delimiter="\t")
    # read sample info from sample sheet
    nsample = samplesheet.shape[0]
    nreads = [0] * nsample
    nunassignedreads = 0
    ntotreads = 0
    SampleIDs = list(samplesheet.SampleID)
    TagReads  = list(samplesheet.TagRead)
    lengthTag = len(TagReads[1])
    print "allowed mismatches = {}".format(options.mismatches)
    
    # check read 1 fastq is present, if so open it
    if not os.path.isfile(options.r1_fn):
        raise IOError("# Error: read 1 file {} does not exist".format(options.r1_fn))
    r1_fgf = gzip.open(options.r1_fn, 'rt')
    
    # check read 2 fastq is supplied and present, if so open it
    is_paired_end = bool(options.r2_fn)
    if (is_paired_end) :
        if not os.path.isfile(options.r2_fn):
            raise IOError("# Error: read 2 file {} does not exist".format(options.r2_fn))
        r2_fgf = gzip.open(options.r2_fn, 'rt')
        
    # make output directory and then file handle for each sample and read end
    try:
        os.mkdir(options.outdir)
    except Exception:
        raise IOError("# Error: output directory {} cannot be created".format(options.outdir))
    nreads_fn = options.outdir + "/nreads.txt"
    if (is_paired_end) :
        # make an output file for each tag and paired read
        r1_out_split_hs = [ gzip.open(options.outdir + "/" + SampleID + "_R1.fastq.gz","wt") \
                for SampleID in SampleIDs]
        r1_out_unassigned_h = gzip.open(options.outdir + "/Unassigned_R1.fastq.gz","wt")
        r2_out_split_hs = [ gzip.open(options.outdir + "/" + SampleID + "_R2.fastq.gz","wt") \
                for SampleID in SampleIDs]
        r2_out_unassigned_h = gzip.open(options.outdir + "/Unassigned_R2.fastq.gz","wt")
    else :
        # make one output file for each tag
        r1_out_split_hs = [ gzip.open(options.outdir + "/" + SampleID + ".fastq.gz","wt") \
            for SampleID in SampleIDs]
        r1_out_unassigned_h = gzip.open(options.outdir + "/Unassigned.fastq.gz","wt")
    
    ## This loop is the heart of the program
    while True:
        # get fastq record/read (4 lines)
        fqrec1 = list(islice(r1_fgf, 4))
        if not fqrec1 :
            break
        if (is_paired_end) :
            fqrec2 = list(islice(r2_fgf, 4))
        
        # count number of processed reads, output every millionth
        ntotreads += 1
        if (ntotreads % 1000000) == 0:
            print "{} reads processed".format(ntotreads)
        
        # assign read to a SampleID
        # TagRead is 1st read with less than threshold mismatches.
        # could cause problems if many mismatches.
        assigned = False
        for sample in range(nsample) :
            # test if initial segment of read matches sample barcode
            if startswith_mismatch(fqrec1[1],TagReads[sample],options.mismatches) :
                assigned = True
                if (options.trimout) :
                    # write trimmed record to file
                    r1_out_split_hs[sample].writelines( trim_fastq_record(fqrec1,n=lengthTag) )
                else : 
                    # write full record to file
                    r1_out_split_hs[sample].writelines( fqrec1 )
                if (is_paired_end) :
                    r2_out_split_hs[sample].writelines( fqrec2 )
                # count it
                nreads[sample] += 1
                # stop testing against other barcodes
                break
        if not assigned: 
            # write unassigned read to file
            # note unassigned reads are not trimmed
            r1_out_unassigned_h.writelines( fqrec1 )
            if (is_paired_end) :
                r2_out_unassigned_h.writelines( fqrec2 )
            # count it
            nunassignedreads += 1
    
    # close output handles and fasta file
    for fh in r1_out_split_hs :
        fh.close()
    r1_out_unassigned_h.close()
    r1_fgf.close()
    if (is_paired_end) :
        for fh in r2_out_split_hs :
            fh.close()
        r2_out_unassigned_h.close()
        r2_fgf.close()
    
    print "All {} reads processed".format(ntotreads)
    
    # output number of reads by sample to file
    samplesheet["nreads"] = nreads
    samplesheet[['SampleID','TagRead','nreads']].to_csv(nreads_fn,sep="\t",index=False)
    # append unassigned reads to samplesheet
    with open(nreads_fn,"a") as nrf :
        nrf.write("Unassd\tNNNNNNNNN\t{}\n".format(nunassignedreads))
        nrf.write("TOTAL\t\t{}".format(ntotreads))
    
    # missing: function call/comments in output
    
    print "Done"


