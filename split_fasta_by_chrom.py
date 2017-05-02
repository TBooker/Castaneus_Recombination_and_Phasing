import sys
from Bio import SeqIO

input_fasta = sys.argv[1]
### Specific to the example of the mm9 build

chrom_tags = []
full_names = []
seq_dict  ={}
with open(input_fasta) as FILE:
	for line in FILE:
		if line.startswith(">"):
			name = line.split(" ")
			if len(name) < 7:
				continue
			if "unlocalized" in name or "unplaced" in name:
				continue
			seq_dict[name[0].strip(">")] = "_".join([name[5],name[6].strip(",")])
			print seq_dict
			#else:
			#	chrom_tags.append(name[0])							
			#	full_names.append(line.strip("\n"))
print "got names list"
fasta = SeqIO.parse(open(input_fasta,"r"),"fasta")
all_chroms = open("all_chromosomes.fasta","w")

for i in fasta:
	if i.id in seq_dict.keys():
		if ">" + i.id in chrom_tags:
			all_chroms.write("\n" + ">"+str(i.id)+"_"+seq_dict[i.id] +"\n" + str(i.seq) + "\n")






		
