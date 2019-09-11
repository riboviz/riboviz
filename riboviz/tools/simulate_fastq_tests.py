#! python

## simulate_fastq_tests.py
## creates simple simulated fastq files to test UMI/deduplication, adapter trimming
## later we can adapt it to test demultiplexing 


# To-do!!
# - Matched functions to make fastq records from a read and quality score plus UMI at 3' end
# -- one places the UMI in-line in the read (raw data, i.e. input for extract; this is done already)
# -- function to places the UMI in the header (i.e. desired output for extract)
# -- these should have the same quality score so we can check quality is preserved by extract
# - Add barcodes for demultiplexing

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from random import choices, seed

QualityMed  = list(range(30,41))
QualityHigh = list(range(39,41))

def simulate_qual(k,QualityVals=QualityMed,weights=None) :
    """
    simulate quality scores; this is a thin wrapper around choices
    whose default values represent medium Phred quality values
    """
    return choices(QualityVals,k=k,weights=weights)

def make_fastq_record(name,readall,qualall=None,QualityVals=QualityMed) :
    """
    make a fastq record with sequence readall and name
    """
    if qualall is None :
        qualall = simulate_qual( len(readall), QualityVals=QualityVals )
    record = SeqRecord(Seq(readall,IUPAC.ambiguous_dna),
               id=name, name=name,
               description=name)
    record.letter_annotations["phred_quality"] = qualall
    return record

if __name__=="__main__" :  
    make_fastq_record("SRR","AAAA")
    seed(42) # fix random seed so test results are the same
    
    ## Components for simulated reads compatible with the vignette files yeast_YAL_CDS_w_250utrs.fa
    ## These are aimed at the Duncan & Mata format with 4nt UMI at each end of read
    ReadA = "ATGGCATCCACCGATTTCTCCAAGATTGAA" # 30nt starting ORF of YAL003W
    ReadAe = "ATGGCATCCACCGATGTCTCCAAGATTGAA" # 1 error in read A
    ReadB = "TCTAGATTAGAAAGATTGACCTCATTAA" # 28nt immediately following start of ORF of YAL038W
    
    UMIX = "CGTA"
    UMIY = "ATAT"
    UMIYe = "ATAA" # 1 error in UMIY
    UMIZ = "CGGC"
    UMIZe = "CTGC" # 1 error in UMIX
    UMI5 = "AAAA"
    UMI5C = "CCCC"
    
    Adaptseq = "CTGTAGGCACC" # adaptor sequence used in vignette data
    
    ### simulate raw data with only 3' UMI, i.e. desired input after adapter trimming
    ### UMI-tools Keeps the highest-quality read, so we use QualityHigh for reads we wish to keep
    sim_records_UMI3_4nt = [ 
    make_fastq_record( "EWSim1ReadAUmiX.Keep",  ReadA + UMIX, QualityVals=QualityHigh), 
    make_fastq_record( "EWSim2ReadAUmiX.Drop",  ReadA + UMIX ), 
    make_fastq_record( "EWSim3ReadAeUmiX.Drop", ReadAe + UMIX ), 
    make_fastq_record( "EWSim4ReadAUmiY.Keep",  ReadA + UMIY ), 
    make_fastq_record( "EWSim5ReadBUmiX.Keep",  ReadB + UMIX ), 
    make_fastq_record( "EWSim6ReadBUmiZ.Keep",  ReadB + UMIZ, QualityVals=QualityHigh), 
    make_fastq_record( "EWSim7ReadBUmiZ.Drop",  ReadB + UMIZ ), 
    make_fastq_record( "EWSim8ReadBUmiZe.Drop", ReadB + UMIZe ) 
    ]
    
    with open("../data/simdata_UMI3_4nt.fastq", "w") as output_handle:
        SeqIO.write(sim_records_UMI3_4nt, output_handle, "fastq")
    
    ### simulate raw data with 3' and 5' UMIs, i.e. desired input after adapter trimming
    sim_records_UMI5and3_4nt = [ 
    make_fastq_record( "EWSim1Umi5ReadAUmiX.Keep",  UMI5 + ReadA + UMIX, QualityVals=QualityHigh), 
    make_fastq_record( "EWSim2Umi5ReadAUmiX.Drop",  UMI5 + ReadA + UMIX ), 
    make_fastq_record( "EWSim3Umi5ReadAeUmiX.Drop", UMI5 + ReadAe + UMIX ), 
    make_fastq_record( "EWSim4Umi5ReadAUmiY.Keep",  UMI5 + ReadA + UMIY ), 
    make_fastq_record( "EWSim5Umi5ReadBUmiX.Keep",  UMI5 + ReadB + UMIX ), 
    make_fastq_record( "EWSim6Umi5ReadBUmiZ.Keep",  UMI5 + ReadB + UMIZ ), 
    make_fastq_record( "EWSim7Umi5ReadBUmiZ.Drop",  UMI5 + ReadB + UMIZ, QualityVals=QualityHigh), 
    make_fastq_record( "EWSim8Umi5ReadBUmiZe.Drop", UMI5 + ReadB + UMIZe ),
    make_fastq_record( "EWSim9Umi5CReadBUmiZ.Keep", UMI5C + ReadB + UMIZe ) 
    ]
    
    with open("../data/simdata_UMI5and3_4nt.fastq", "w") as output_handle:
        SeqIO.write(sim_records_UMI5and3_4nt, output_handle, "fastq")
    
    ### simulate raw data with 3' and 5' UMIs, i.e. desired input after adapter trimming
    sim_records_UMI5and3_4nt_adaptor = [ 
    make_fastq_record( "EWSim1Umi5ReadAUmiX.Keep",  UMI5 + ReadA + UMIX + Adaptseq, QualityVals=QualityHigh ), 
    make_fastq_record( "EWSim2Umi5ReadAUmiX.Drop",  UMI5 + ReadA + UMIX + Adaptseq ), 
    make_fastq_record( "EWSim3Umi5ReadAeUmiX.Drop", UMI5 + ReadAe + UMIX + Adaptseq ), 
    make_fastq_record( "EWSim4Umi5ReadAUmiY.Keep",  UMI5 + ReadA + UMIY + Adaptseq ), 
    make_fastq_record( "EWSim5Umi5ReadBUmiX.Keep",  UMI5 + ReadB + UMIX + Adaptseq + "AC"), # extra nt past the adaptor for the shorter read
    make_fastq_record( "EWSim6Umi5ReadBUmiZ.Keep",  UMI5 + ReadB + UMIZ + Adaptseq + "AC", QualityVals=QualityHigh), 
    make_fastq_record( "EWSim7Umi5ReadBUmiZ.Drop",  UMI5 + ReadB + UMIZ + Adaptseq + "AC" ), 
    make_fastq_record( "EWSim8Umi5ReadBUmiZe.Drop", UMI5 + ReadB + UMIZe + Adaptseq + "AC" ),
    make_fastq_record( "EWSim9Umi5CReadBUmiX.Keep", UMI5C + ReadB + UMIX + Adaptseq + "AC" ) 
    ]
    
    with open("../data/simdata_UMI5and3_4nt_adaptor.fastq", "w") as output_handle:
        SeqIO.write(sim_records_UMI5and3_4nt_adaptor, output_handle, "fastq")
