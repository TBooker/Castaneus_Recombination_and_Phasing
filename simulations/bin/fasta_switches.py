import sys, random
from Bio import SeqIO
## This script takes a fasta file of haplotypes, assumes that an individuals chromosomes are named i, i+1
##
## It takes a switch error rate to assess whether simulations 
##
def wrap(to_wrap, width=70):
    return '\n'.join(to_wrap[i:i+width] for i in range(0, len(to_wrap), width))

if len(sys.argv) < 4:
	print "need to ive the following commands at the command line:\nPython fasta_switches.py INPUT.fa SWITCH_RATE OUTPUT.fa\n"
	sys.exit()

fasta = SeqIO.parse(open(sys.argv[1],"r"),"fasta")
switch_rate = float(sys.argv[2])
output = open(sys.argv[3],"w")


seq_dict = {}
for i in fasta:
#	print i.id

	
	uniq = i.id.split("_")[0]
	id_num = int(i.id.split("_")[-2])
	name =  uniq+"_"+str(id_num)
	#continue
	seq_dict[name] = str(i.seq)

individuals = []

for i in seq_dict.keys():
#	print i
	uniq = i.split("_")[0]
	id_num = int(i.split("_")[-1])
	if id_num % 2 == 0:
		individuals.append([uniq+"_"+str(id_num-1),uniq+"_"+str(id_num)])
#print
#print seq_dict.keys()
#print
#max_length =0
for ind in individuals:
#	print ind
#	print ind[0],len(seq_dict[ind[0]])
#	print ind[1],len(seq_dict[ind[1]])

	switch = 1
	het_count = 0
	switch_count =0
	_1_ = "" #the strings for the new individuals
	_2_ = ""
	for i,j in zip(seq_dict[ind[0]],seq_dict[ind[1]]):
#		max_length +=1
		if i != j:
			het_count += 1
			switch_prob =random.random()
			if switch_prob <= switch_rate:
				switch_count+=1
				switch *= -1
		if switch < 0:
			_1_ += j
			_2_ += i
		if switch > 0:
			_1_ += i
			_2_ += j
#		if max_length == 50:
#			max_length = 0 
#			break


	print ind, het_count, switch_count
	output.write(">"+ind[0]+"\n"+wrap(_1_)+"\n")
	output.write(">"+ind[1]+"\n"+wrap(_2_)+"\n")

output.close()
