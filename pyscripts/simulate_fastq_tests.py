#! python

## simulate_fastq_tests.py
## creates simple simulated fastq files to test UMI/deduplication, adapter trimming
## later we can adapt it to test demultiplexing 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from random import choices

ReadA = "ATGGCATCCACCGATTTCTCCAAGATTGAA" # 30nt starting ORF of YAL003W
ReadB = "TCTAGATTAGAAAGATTGACCTCATTAA" # 28nt immediately following start of ORF of YAL038W

UMIX = "CGTA"
UMIY = "ATAT"
UMIYe = "ATAA"
UMIZ = "CGGC"
UMIZe = "CTGC"

# we will use QualityAlphabet for random sampling of phred scores
QualityAlphabet = [ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J" ]

def simulate_qual(k) :
	return "".join( choices(QualityAlphabet,k=k) )

def make_fastq_record(name,readall,qualall=None) :
	"""
	make a fastq record with sequence readall and name
	"""
	if qualall is None :
		qualall = simulate_qual( len(readall) )
	record = SeqRecord(Seq(readall,IUPAC.ambiguous_dna),
			   id=name, name=name,
			   description=name)
	return record

make_fastq_record("SRR","AAAA")

# To-do!!
# - Check that make_fastq_record produces output writeable into a fastq file by SeqIO.write()
# - Add quality scores to fastq record
# - Two functions to make fastq records from a read and quality score plus UMI at 3' end
# -- one places the UMI in-line in the read (raw data, i.e. input for extract)
# -- the other places the UMI in the header (i.e. desired output for extract)
# -- these should have the same quality score so we can check quality is preserved by extract
# - Then we can produce a desired input and output fastq file with ReadA and ReadB (on board, 27th August)
# - Then extend to UMIs at BOTH 5' and 3' end of read
# - After that we can test adapter trimming and/or barcodes for demultiplexing