#! python

## Check fasta and gff files are compatible with RiboViz
## Specifically, test that 
##    - the beginning of every CDS is a start codon (ATG; translates to M)
##    - the stop of every CDS is a stop codon (TAG, TGA, TAA; translates to *)
##    - there are no stop codons internal to the CDS.
## Some unusual genes (e.g. frameshifts) might not have this.
## 
## example:
##   python scripts/check_fasta_gff.py -fa vignette/input/yeast_YAL_CDS_w_250utrs.fa -gff vignette/input/yeast_YAL_CDS_w_250utrs.gff3 

import argparse, gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from warnings import warn

if __name__=="__main__" :
    # take input options
    parser = argparse.ArgumentParser(description="Check fasta and gff files have start and stop codons as expected")
    parser.add_argument("-fa","--fastain",dest="fastain",nargs='?',help="fasta file input")
    parser.add_argument("-gff","--gffin",dest="gffin",nargs='?',help="gff3 file input")
    options = parser.parse_args()
    
    gffin   = options.gffin
    fastain = options.fastain
    
    print( "Checking fasta file " + fastain + "\n with gff file " + gffin )
    
    # open fasta and gff files
    gffdb = gffutils.create_db(gffin, dbfn='test.db', force=True, keep_order=True,
    merge_strategy='merge', sort_attribute_values=True)
    
    for CDS_coord in gffdb.features_of_type('CDS'):
        # for all CDS entries in gff,
        # print(CDS_coord.seqid)
        # extract CDS
        CDS_seq = CDS_coord.sequence(fastain)
        
        CDS_len_remainder = len(CDS_seq) % 3 
        if not ( CDS_len_remainder == 0 ) :
           warn( CDS_coord.seqid + " has length that isn't divisible by 3" )
           CDS_seq += ( "N" * ( 3 - CDS_len_remainder) )
        
        # translate CDS
        CDS_trans = Seq( CDS_seq, IUPAC.ambiguous_dna ).translate()
        
        if ( CDS_trans[0] != "M" ) :
            warn( CDS_coord.seqid + " doesn't start with ATG" )
        if ( CDS_trans[-1] != "*" ) :
            warn( CDS_coord.seqid + " doesn't stop at end" )
        if any( [ L == "*"  for L in CDS_trans[:-1] ] ) : 
            warn( CDS_coord.seqid + " has internal STOP" )
    
    print("Done checking")