#!/usr/bin/env python3

#gffutils, biopython, bcbio-gff packages need to be installed
import gffutils,os,argparse
from Bio import SeqIO
from collections import defaultdict
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


parser = argparse.ArgumentParser()
parser.add_argument('--fasta')
parser.add_argument('--gff3')
parser.add_argument('--output')
parser.add_argument('--min_length')
args = parser.parse_args()

myfasta = args.fasta
mygff3 = args.gff3
output_file = args.output
protein_length = int(args.min_length)

def getORF(seq, treshold, start, stop, end_UTR5, end_CDS, name):
    f = defaultdict(list)
    for frame in range(3):
         start_codon_index = 0
         end_codon_index = 0
         start_codon_found = False
         for i in range(frame, len(seq), 3):
              current_codon = seq[i:i+3]                            
              if current_codon in start and not start_codon_found:
                 start_codon_found = True
                 start_codon_index = i
                 e = start_codon_index + 1                
              if current_codon in stop and start_codon_found:
                 end_codon_index = i
                 d = end_codon_index+3
                 length = end_codon_index - start_codon_index + 3
                 if length >= treshold * 3:        #treshold is protein length but length is nuc
                     b = length//3 + 3
                     start_codon_found = False
                     if e < end_UTR5 and d < end_UTR5:
                        type = "uORF"
                        list1 = [b, length, frame+1, e, d, type, name]
                        f[name].append(list1)
                        #print("aa %i, nuc %i, frame %i, coord %i:%i type %s name %s\n" % (b, length, frame+1, e, d, type, name))
                        #f.write(name + '\t' + type + '\t' + str(e) +'\t'+ str(d)+ '\n')
                     if e < end_UTR5 and d == end_CDS:
                        type = "CDS_NTE"
                        list2 = [b, length, frame+1, e, d, type, name]
                        f[name].append(list2)
                        #print("aa %i, nuc %i, frame %i, coord %i:%i type %s name %s\n" % (b, length, frame+1, e, d, type, name))
                        #f.write(name + '\t' + type + '\t' + str(e) + '\t'+ str(d)+ '\n')
                     if e < end_UTR5 and end_UTR5 < d < end_CDS:
                        type = "overlap_uORF"
                        list2 = [b, length, frame+1, e, d, type, name]
                        f[name].append(list2)
                        #print("aa %i, nuc %i, frame %i, coord %i:%i type %s name %s\n" % (b, length, frame+1, e, d, type, name))   
                        #f.write(name + '\t' + type + '\t' + str(e) +'\t'+ str(d)+ '\n')
    return f
    
db = gffutils.create_db(mygff3, 'myGFF.db', merge_strategy="create_unique", keep_order=True, force=True)
db = gffutils.FeatureDB('myGFF.db')

position_dict = {}
for i in db.features_of_type("UTR5"):
    dict = {}
    for j in db.features_of_type("CDS"):
         if i.seqid == j.seqid:
            dict[i.seqid]=[i.start-1,j.end,i.end]
            position_dict.update(dict)
            
#start codon
start = ["ATG","GTG","TTG"]

#stop codon
stop = ["TAA","TAG","TGA"]

for key in position_dict.keys():
    
    for record in SeqIO.parse(myfasta, "fasta"):
        
               if key == record.id:
                 
                  sequence = record.seq[position_dict[key][0]:position_dict[key][1]]
                 
                  end_UTR5 = position_dict[key][2]
                  
                  end_CDS = position_dict[key][1]
                 
                  name = key
                  
                  print(name)
                  
                  treshold = protein_length
                 
                  output_dict = getORF(sequence,treshold,start,stop,end_UTR5,end_CDS,name)
                  
                  print(output_dict)
                
                  for values in output_dict.values():
                        
                        print(values[0][3])
                        
                        for i in range(len(values)):
                            
                             out_file = output_file
                  
                             seq = record.seq
                  
                             rec = SeqRecord(seq, name)
                  
                             qualifiers = {"source": "riboviz", "score": ".", "start_codon": seq[values[i][3]-1:values[i][3]+2], "Name": name + "_" + values[i][5] + "_" + str(values[i][3]),"frame":values[i][2]}
                  
                             feature = SeqFeature(FeatureLocation(values[i][3]-1, values[i][4]), type=values[i][5], strand=1,qualifiers=qualifiers)

                             rec.features = [feature]

                             with open(out_file, "a") as out_handle:
                      
                                GFF.write([rec], out_handle)
                                
os.system("sed -i '/#/d' {}" .format(output_file))
