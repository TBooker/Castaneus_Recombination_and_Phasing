### This script will generate random ATCG sequence for a given SLiM run and convert the slim run into a fasta file for a random individual in the population

import random, argparse,tom_slim, glob, collections
import cPickle as pickle
from tom import brace

def mutate(original,mut_freqs):	
	alt = random.choice("ATCG")
	random_position = random.uniform(0,1)
	proportions = [mut_freqs[original+":"+j] for j in "ACGT"]
	accumulator = 0
	for i,j in zip(proportions,"ACGT"):
		accumulator += i
		if random_position > accumulator - i and random_position < accumulator :
			base_selected = j
	return base_selected


def get_ref(mut_mat):
	lengths = [mut_mat[i+":."] for i in "ACGT"]
	proportions = [k/sum(lengths) for k in lengths]
	random_position = random.uniform(0,1)
	accumulator = 0
	for i,j in zip(proportions,"ACGT"):
		accumulator += i
		if random_position > accumulator - i and random_position < accumulator :
			base_selected = j
	return base_selected	





def generate_seq_dicts(slim_in,mat_dict,number_of_sequences,all_individuals=True):
	slim = tom_slim.slim(slim_in,fixed=False,give_genomes=True,all_individuals=True)
	mut_mat = mat_dict["mut_mat"]
	mut_freqs = mat_dict["mut_freqs"]
	reference =  tuple([get_ref(mut_mat) for i in xrange(slim.length)]) ## This defines the ancestral sequence as a string of random letters (ATCG) with the length of the SLiM chromosome
#	return
	#### LEFT IT HERE JANUARY 20th. Need to make a dict of the mutation frequencies so that the mutate function can do the same thing as the get_ref function, to choose the base to mutate to based upon the matrix

	mut_dict_raw = slim.mutations_dict(minFreq=1)
	
#	for i in mut_dict_raw:
#		print i, mut_dict_raw[i]
	mut_dict = {}

	for i in mut_dict_raw.keys():
		pos = mut_dict_raw[i]-1
		ref_at_pos = reference[pos]
		alt_at_pos = mutate(ref_at_pos,mut_freqs)
		mut_dict[i] = [pos,ref_at_pos,alt_at_pos]

	samples = slim.genome_dict()
	
	#haps = {}
	haps = collections.OrderedDict()
	individuals_chosen = []
	while len(individuals_chosen) < number_of_sequences:
	### Add a loop here to get multiple individuals
		individual_to_choose = random.randint(1,(len(samples.keys())/2))
		if individual_to_choose in individuals_chosen: continue
		else: individuals_chosen.append(individual_to_choose)
		seq_1 = "p1:"+str(individual_to_choose *2)
		seq_2 = "p1:"+str(individual_to_choose *2 -1)
#		brace()
		
		for p in [seq_1,seq_2]:
			haps[p] = samples[p]
		
	return haps,slim.name,mut_dict

def wrap(to_wrap, width=70):
    return '\n'.join(to_wrap[i:i+width] for i in range(0, len(to_wrap), width))

def random_sample(total,num_samples):
	sample = []
	while len(sample) != num_samples:
		x = random.randint(1,total)
		if x not in sample:
			sample.append(x)
	return sample




parser = argparse.ArgumentParser(description="This script iterates through a slim output file and generates a dictionary which can be used at a later point to generate a set of FASTA sequences for further analyses")

parser.add_argument("-i","--input", 
		required = True,
		dest = "input", 
		type =str, 
		help = "The SLiM output file")
parser.add_argument("-o","--output", 
		required = False,
		dest = "output", 
		type =str, 
		help = "The name of the output dictionary")
parser.add_argument("--gz", 
		required = False,
		dest = "gz", 
		action = "store_true", 
		help = "Is the input file gzipped? It is recommended (by me) that the file should be gzipped",
		default = False)
parser.add_argument("--procs","-p", 
		required = False,
		dest = "procs", 
		type = int, 
		help = "How many cores do you want to engage in this operation?",
		default = 1)
parser.add_argument("--sequences","-n", 
		required = False,
		dest = "sequences", 
		type = int, 
		help = "How many individuals' sequences do you to extract?",
		default = 1)
parser.add_argument("--mut_matrix","-m", 
		required = True,
		dest = "mut_mat", 
		type = str, 
		help = "Give the name of the mutation matrix that you want to use")


args = parser.parse_args()


## Get the default output name if not supplied
if not args.output:
	if args.gz:
		output_name = args.input.replace(".gz",".seq.pkl")
	else:
		output_name = args.input +".seq.pkl"
else:
	output_name = args.output+".seq.pkl"


## This is a new idea for me, if the user (me in the future, wants to extact the SFS for a gzipped set of slim runs, they can.  This little bed gets the right function for the job
if args.gz:
	slim_iterator = tom_slim.slim_reader_gzip
else:
	slim_iterator = tom_slim.slim_reader
output_dict = {}
# put all the dictionaries here so that they are only being defined once
mut_mat = {}
for k in open(args.mut_mat,"r"):
	cc = k.strip("\n").split(" ")
	mut_mat[cc[0]] = float(cc[1])


mut_mat_keys = [i for i in mut_mat.keys() if i[2] != "."]
mut_freqs = {}
mut_counts = {}
for l in "ACGT":
	mut_counts[l] = sum([mut_mat[l+":"+p] for p in  "ACGT"])
for l in "ACGT":
	for k in "ACGT":
		mut_freqs[l+":"+k] = mut_mat[l+":"+k]/mut_counts[l]

mat_dict ={"mut_mat":mut_mat,"mut_freqs":mut_freqs}


print "Iterating over all SLiM runs"

haps = {}
for i in slim_iterator(args.input):
	haps_raw,name,mut_dict = generate_seq_dicts(i,mat_dict,args.sequences)
	sample_muts =  set([item for inner_list in [haps_raw[i] for i in haps_raw.keys()] for item in inner_list])
	print sample_muts
	print haps_raw
