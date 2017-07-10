#! python

## Remove a single 5' mismatched nt AND filter reads with more than specified mismatches from sam file
## only checks positive-strand alignments
## example:
##   python trim_5p_mismatch.py -in testdata_trim_5p_mismatch.sam -out test_trim_5p_mismatch_clean.sam
##   python trim_5p_mismatch.py -in data_map1.sam -out data_map1_clean.sam

import pysam, argparse, re

if __name__=="__main__" :  
    parser = argparse.ArgumentParser(description="Remove a single 5' mismatched nt AND filter reads with more than specified mismatches from sam file")
    parser.add_argument("-in","--samfilein",dest="samfilein",nargs='?',help="sam file input")
    parser.add_argument("-out","--samfileout",dest="samfileout",nargs='?',help="sam file output")
    parser.add_argument("-mm","--mismatches",dest="mismatches",nargs='?',default=1L,type=int,help="number of mismatches to allow")
    parser.add_argument("-5p","--fivepremove",dest="fivepremove",nargs='?',default=True,help="remove 5' mismatched nt")
    options = parser.parse_args()
    
    print "trim_5p_mismatch.py running"
    
    # open input and output files
    samin=pysam.AlignmentFile(options.samfilein, "r")
    samout=pysam.AlignmentFile(options.samfileout,"wh",template=samin)
    
    # count number of reads processed, discarded, trimmed, written
    nprocessed = 0
    ndiscarded = 0
    ntrimmed   = 0
    nwritten   = 0
    
    for read in samin.fetch():
        # loop over input reads
        nprocessed += 1
        if (nprocessed % 1000000 ) == 1 :
            print "processed " + str(nprocessed - 1) + " reads"
        try :
            # get MD tag for read, encoding mismatches
            MDtag=read.get_tag('MD')
        except KeyError :
            # MD tag not present, assume read not aligned, discard
            ndiscarded += 1
            continue
        # count mismatches in read
        nmismatch = sum( [ MDtag.count(L) for L in "ATCG" ] )
        if options.fivepremove & ( MDtag[0]=="0" ) :
            # If the 5' nt is mismatched... 
            if MDtag[2] in ["A","T","C","G","0"] :
                # 2nd nt is also mismatched; discard
                ndiscarded += 1
                continue
            # ... soft-clip 5' nt
            # increment position of alignment
            read.pos += 1
            # edit MD tag to remove leading mismatch
            read.set_tag('MD', re.sub( "^0[ATCG]", "", MDtag ) )
            # lower the number of mismatches by 1
            nmismatch -= 1
            read.set_tag('NM', nmismatch)
            # edit CIGAR string to increase soft clip
            cigarstring = read.cigarstring
            if ( not "S" in cigarstring ) or ( cigarstring.find("S") > cigarstring.find("M") ) :
                # read is not soft-clipped on left
                # find number of initial matches
                ninitmatch = int( re.findall( r"^([0-9]+)M", cigarstring )[0] )
                # add initial "1S" and reduce initial match by 1
                newinitbit = "1S" + str( ninitmatch - 1L ) + "M"
                read.cigarstring = re.sub( r"^([0-9]+)M", newinitbit, cigarstring )
                ntrimmed += 1
            elif ( cigarstring.find("S") < cigarstring.find("M") ) :
                # read is soft-clipped on left
                nsoftclip = int( re.findall( r"^([0-9]+)S", cigarstring )[0] )
                ninitmatch = int( re.findall( r"([0-9]+)M", cigarstring )[0] )
                # add 1 to left soft-clip and reduce initial match by 1
                newinitbit = str( nsoftclip + 1L ) + "S" + str( ninitmatch - 1L ) + "M"
                read.cigarstring = re.sub( r"^([0-9]+)S([0-9]+)M", newinitbit, cigarstring )
                ntrimmed += 1
            else :
                # Cry for help and discard
                print "Odd cigar string for read : "
                print read
                ndiscarded += 1
                continue
        if nmismatch <= options.mismatches : 
            # write to output if <= mm mismatches
            nwritten += 1
            samout.write(read)
        else :
            ndiscarded += 1
            
    # close input and output files
    samin.close()
    samout.close()
    
    # report read counts and exit
    print "trim_5p_mismatch.py finished, number of reads"
    print "processed:\t" + str(nprocessed)
    print "discarded:\t" + str(ndiscarded)
    print "trimmed:\t" + str(ntrimmed)
    print "written:\t" + str(nwritten)
