## a script to get the mutations that can be used to infer the mutation rate matrix used by LDhelmet

import gzip, argparse

parser = argparse.ArgumentParser(description="Takes an annotation table and a .Freq file and generates a list of sites of a particular fold")

parser.add_argument("-c", "--chr", 
		required = True, 
		metavar ="chr", 
		type =str, 
		help = "The chromosome that you want to use (Don't prefix with chr)")
parser.add_argument("--cpg", 
		required = False, 
		metavar ="cpg", 
		type=int, 
		help = "Do you want CpGs? Give 0 for no, 1 for both CpG and non and 2 for only CpG sites",
		default = 0)

args = parser.parse_args()


freq_file = "/home/booker/mouse_genome/freq_files/chr"+args.chr+"/chr"+args.chr+".freq.gz"

seq_dict = {0:"A",1:"C",2:"G",3:"T"}

mut_dict = {}

base_dict = {}

for i in "ACGT":
	base_dict[i] = 0
	mut_dict[i] = {}
	for j in "ACGT":
		mut_dict[i][j] = 0
#print mut_dict
#count = 0
for p in gzip.open(freq_file):
#	if count == 10000:break
	line = p.strip("\n").split("\t")
	
	CpG_status = set([i for i in line[4],line[6],line[8] if i != "."])
	if args.cpg ==1: 
		pass
	else:
		if args.cpg == 0 and CpG_status !=  set(["0"]): 
				continue
	
		elif args.cpg == 2 and CpG_status !=  set(["1"]): 
				continue
	try:
		fam_allele = seq_dict[line[7].split(",").index("2")]
	except ValueError:continue
	try:
		rat_allele = seq_dict[line[9].split(",").index("1")]
	except ValueError:continue

	cast = line[5].split(",")
	cast_alleles = [seq_dict[i] for i in xrange(len(cast)) if cast[i]!="0"]
	
	if len(cast_alleles) == 1:
		base_dict[cast_alleles[0]] +=1
		continue
	elif len(cast_alleles) == 3: continue
	else: pass
		 ## gets bi-allelic sites

	try:
		fam_allele = seq_dict[line[7].split(",").index("2")]
	except ValueError:continue
	
	rat_allele = seq_dict[line[9].split(",").index("1")]

	if rat_allele != fam_allele:continue ## Ancestral allele is now rat_allele
	
	alt_allele = [i for i in cast_alleles if i != rat_allele][0]
	mut_dict[rat_allele][alt_allele] += 1
#	count +=1

output = open(args.chr+"_mut_rates.txt","w")
for i in mut_dict.keys():
	output.write(i+"\t"+"."+"\t"+str(base_dict[i])+"\n")
	for j in mut_dict[i].keys():
		output.write(i+"\t"+j+"\t"+str(mut_dict[i][j])+"\n")

output.close()	
