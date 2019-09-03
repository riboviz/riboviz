#! python

## simulate_fastq_tests.py
## creates simple simulated fastq files to test UMI/deduplication, adapter trimming
## later we can adapt it to test demultiplexing 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from random import choices

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

if __name__=="__main__" :  
	make_fastq_record("SRR","AAAA")
	
	#### Components for simulated reads compatible with the vignette files yeast_YAL_CDS_w_250utrs.fa
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
	
	sim_records_UMI3_4nt <- [ 
	make_fastq_record( "EW_Sim1_ReadA_UX_Keep",  ReadA + UMIX ), 
	make_fastq_record( "EW_Sim2_ReadA_UX_Drop",  ReadA + UMIX ), 
	make_fastq_record( "EW_Sim3_ReadAe_UX_Drop", ReadAe + UMIX ), 
	make_fastq_record( "EW_Sim4_ReadA_UY_Keep",  ReadA + UMIY ), 
	make_fastq_record( "EW_Sim5_ReadB_UX_Keep",  ReadB + UMIX ), 
	make_fastq_record( "EW_Sim6_ReadB_UZ_Keep",  ReadB + UMIZ ), 
	make_fastq_record( "EW_Sim7_ReadB_UZ_Drop",  ReadB + UMIZ ), 
	make_fastq_record( "EW_Sim8_ReadB_UZe_Drop", ReadB + UMIZe ) 
	]
	
	with open("../data/simdata_UMI3_4nt.fastq", "w") as output_handle:
    	SeqIO.write(sim_records_UMI3_4nt, output_handle, "fastq")
	
	
	sim_records_UMI5and3_4nt <- [ 
	make_fastq_record( "EW_Sim1_U5_ReadA_UX_Keep",  UMI5 + ReadA + UMIX ), 
	make_fastq_record( "EW_Sim2_U5_ReadA_UX_Drop",  UMI5 + ReadA + UMIX ), 
	make_fastq_record( "EW_Sim3_U5_ReadAe_UX_Drop", UMI5 + ReadAe + UMIX ), 
	make_fastq_record( "EW_Sim4_U5_ReadA_UY_Keep",  UMI5 + ReadA + UMIY ), 
	make_fastq_record( "EW_Sim5_U5_ReadB_UX_Keep",  UMI5 + ReadB + UMIX ), 
	make_fastq_record( "EW_Sim6_U5_ReadB_UZ_Keep",  UMI5 + ReadB + UMIZ ), 
	make_fastq_record( "EW_Sim7_U5_ReadB_UZ_Drop",  UMI5 + ReadB + UMIZ ), 
	make_fastq_record( "EW_Sim8_U5_ReadB_UZe_Drop", UMI5 + ReadB + UMIZe ),
	make_fastq_record( "EW_Sim9_U5C_ReadB_UZ_Keep", UMI5C + ReadB + UMIZe ) 
	]
	
	with open("../data/simdata_UMI5and3_4nt.fastq", "w") as output_handle:
    	SeqIO.write(sim_records_UMI5and3_4nt, output_handle, "fastq")
	
	sim_records_UMI5and3_4nt_adaptor <- [ 
	make_fastq_record( "EW_Sim1_U5_ReadA_UX_Keep",  UMI5 + ReadA + UMIX + Adaptseq ), 
	make_fastq_record( "EW_Sim2_U5_ReadA_UX_Drop",  UMI5 + ReadA + UMIX + Adaptseq ), 
	make_fastq_record( "EW_Sim3_U5_ReadAe_UX_Drop", UMI5 + ReadAe + UMIX + Adaptseq ), 
	make_fastq_record( "EW_Sim4_U5_ReadA_UY_Keep",  UMI5 + ReadA + UMIY + Adaptseq ), 
	make_fastq_record( "EW_Sim5_U5_ReadB_UX_Keep",  UMI5 + ReadB + UMIX + Adaptseq + "AC"), # extra nt past the adaptor for the shorter read
	make_fastq_record( "EW_Sim6_U5_ReadB_UZ_Keep",  UMI5 + ReadB + UMIZ + Adaptseq + "AC" ), 
	make_fastq_record( "EW_Sim7_U5_ReadB_UZ_Drop",  UMI5 + ReadB + UMIZ + Adaptseq + "AC" ), 
	make_fastq_record( "EW_Sim8_U5_ReadB_UZe_Drop", UMI5 + ReadB + UMIZe + Adaptseq + "AC" ),
	make_fastq_record( "EW_Sim9_U5C_ReadB_UX_Keep", UMI5C + ReadB + UMIX + Adaptseq + "AC" ) 
	]
	
	with open("../data/simdata_UMI5and3_4nt_adaptor.fastq", "w") as output_handle:
    	SeqIO.write(sim_records_UMI5and3_4nt_adaptor, output_handle, "fastq")
	


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