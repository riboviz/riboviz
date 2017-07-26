#! python

## Remove a single 5' mismatched nt AND filter reads with more than specified mismatches from sam file
## example:
##   python trim_5p_mismatch.py -in testdata_trim_5p_mismatch.sam -out testdata_trim_5p_mismatch_clean.sam
##   python trim_5p_mismatch.py -in testdata_trim_5pos5neg.sam -out testdata_trim_5pos5neg_clean.sam
##   python trim_5p_mismatch.py -in data_map1.sam -out data_map1_clean.sam

import pysam, argparse, re

if __name__=="__main__" :  
    parser = argparse.ArgumentParser(description="Remove a single 5' mismatched nt AND filter reads with more than specified mismatches from sam file")
    parser.add_argument("-in","--samfilein",dest="samfilein",nargs='?',help="sam file input")
    parser.add_argument("-out","--samfileout",dest="samfileout",nargs='?',help="sam file output")
    parser.add_argument("-mm","--mismatches",dest="mismatches",nargs='?',default=1L,type=int,help="number of mismatches to allow")
    fivep_parser = parser.add_mutually_exclusive_group(required=False)
    fivep_parser.add_argument("-5p","--fivepremove", dest='fivepremove', action='store_true')
    fivep_parser.add_argument("-no5p","--nofivepremove", dest='fivepremove', action='store_false')
    parser.set_defaults(fivepremove=True)
    options = parser.parse_args()
    print options
    
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
        nmismatch = read.get_tag('NM')
        # nmismatch = sum( [ MDtag.count(L) for L in "ATCG" ] )
        
        if nmismatch > 0 & options.fivepremove:
            # if there are any mismatches ..
            
            # import pdb; pdb.set_trace()
            if ( MDtag[0]=="0" ) & (read.flag==0) :
                # If the 5' nt is mismatched on a plus-strand read... 
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
                if ( not "S" in cigarstring[:3] ) :
                    # read is not soft-clipped on left
                    # find number of initial matches
                    ninitmatch = int( re.findall( r"^([0-9]+)M", cigarstring )[0] )
                    # add initial "1S" and reduce initial match by 1
                    newinitbit = "1S" + str( ninitmatch - 1L ) + "M"
                    read.cigarstring = re.sub( r"^([0-9]+)M", newinitbit, cigarstring )
                else :
                    # read is soft-clipped on left
                    nsoftclip = int( re.findall( r"^([0-9]+)S", cigarstring )[0] )
                    ninitmatch = int( re.findall( r"([0-9]+)M", cigarstring )[0] )
                    # add 1 to left soft-clip and reduce initial match by 1
                    newinitbit = str( nsoftclip + 1L ) + "S" + str( ninitmatch - 1L ) + "M"
                    read.cigarstring = re.sub( r"^([0-9]+)S([0-9]+)M", newinitbit, cigarstring )
                # count the read as trimmed
                ntrimmed += 1
        
            if ( bool( re.search( "[ATCG]0$", MDtag ) ) ) & (read.flag==16) :
                # If the 5' nt is mismatched on a minus strand read... 
                # positive sense is with template; read is reverse-complement
                if MDtag[-3] in ["A","T","C","G"] :
                    # 2nd nt is also mismatched; discard
                    ndiscarded += 1
                    continue
                # ... soft-clip 5' nt
                # don't increment position of alignment!
                # edit MD tag to remove trailing mismatch
                read.set_tag('MD', re.sub( "[ATCG]0$", "", MDtag ) )
                # lower the number of mismatches by 1
                nmismatch -= 1
                read.set_tag('NM', nmismatch)
                # edit CIGAR string to increase soft clip
                cigarstring = read.cigarstring
                if ( cigarstring[-1] != "S" ) :
                    # read is not soft-clipped on left ( = right along template)
                    # find number of terminal matches
                    ntermmatch = int( re.findall( r"([0-9]+)M$", cigarstring )[0] )
                    # add initial "1S" and reduce initial match by 1
                    newtermbit = str( ntermmatch - 1L ) + "M1S"
                    read.cigarstring = re.sub( r"([0-9]+)M$", newtermbit, cigarstring )
                else :
                    # read is soft-clipped on left ( = right along template)
                    # find number of terminal matches
                    nsoftclip = int( re.findall( r"([0-9]+)S$", cigarstring )[0] )
                    ntermmatch = int( re.findall( r"([0-9]+)M", cigarstring )[-1] )
                    # add initial "1S" and reduce initial match by 1
                    newtermbit = str( ntermmatch - 1L ) + "M" + str( nsoftclip + 1L ) + "S"
                    read.cigarstring = re.sub( r"([0-9]+)M([0-9]+)S$", newtermbit, cigarstring )
                # count the read as trimmed
                ntrimmed += 1
        
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
