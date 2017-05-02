import sys, random
from Bio import SeqIO
## This script takes a fasta file of haplotypes, assumes that an individuals chromosomes are named i, i+1
##
## It takes a switch error rate to assess whether simulations 
##
def wrap(to_wrap, width=70):
    return '\n'.join(to_wrap[i:i+width] for i in range(0, len(to_wrap), width))

if len(sys.argv) < 3:
	print "need to ive the following commands at the command line:\nPython fasta_to_sequenceLDhot.py INPUT.fa OUTPUT.file\n"
	sys.exit()

seq_dict = SeqIO.to_dict(SeqIO.parse(open(sys.argv[1],"r"),"fasta"))

output = open(sys.argv[2],"w")

haps = {}
haps["pos"] = []
count = 0
print "getting dict of haplotypes"
for x in range(len(seq_dict[seq_dict.keys()[0]])):

	sites = set([seq_dict[ind].seq[x] for ind in seq_dict.keys()]) 
	if len(sites) > 1:
		count +=1
		haps["pos"].append(str(x + 1))

		for idv in seq_dict.keys():
			if idv not in haps.keys():
				haps[idv] = [seq_dict[idv].seq[x]]
			else:
				haps[idv].append(seq_dict[idv].seq[x])
pos_string =  " ".join(haps["pos"])
print "finding distinct sequences"
seqs = ["".join(haps[i]) for i in haps.keys() if i!="pos"] 
distinct = set(seqs)
print "writing haps to file"
y = str("Distinct ="+ str(len(distinct)) + \
 "\nGenes = "+str(len(seqs))+ \
"\nLoci = "+str(count) + \
"\nK = -4" + \
"\nPositions of loci:\n" + \
pos_string + \
"\nHaplotypes\n")
output.write(y)
[output.write("\t"+sequence +"\t1\n") for sequence in distinct]
		
output.write("#")
output.close()
