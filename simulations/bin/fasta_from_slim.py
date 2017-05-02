### This script will generate random ATCG sequence for a given SLiM run and convert the slim run into a fasta file for a random individual in the population

import random, argparse,tom_slim, glob
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





def generate_seq_dicts(slim_in,mat_dict,number_of_sequences,all_individuals=True,singletons=False):
	slim = tom_slim.slim(slim_in,fixed=True,give_genomes=True,all_individuals=True)

	mut_mat = mat_dict["mut_mat"]
	mut_freqs = mat_dict["mut_freqs"]
	reference =  tuple([get_ref(mut_mat) for i in xrange(slim.length)]) ## This defines the ancestral sequence as a string of random letters (ATCG) with the length of the SLiM chromosome
#	return
	print slim.name
	#### LEFT IT HERE JANUARY 20th. Need to make a dict of the mutation frequencies so that the mutate function can do the same thing as the get_ref function, to choose the base to mutate to based upon the matrix
        if singletons==True:
                min_frequency=0
        elif singletons==False:
                min_frequency=1

	mut_dict_raw = slim.mutations_dict(minFreq=min_frequency)
	
	for i in mut_dict_raw:
		print i, mut_dict_raw[i]
	mut_dict = {}

	for i in mut_dict_raw.keys():
		pos = mut_dict_raw[i]-1
		ref_at_pos = reference[pos]
		alt_at_pos = mutate(ref_at_pos,mut_freqs)
		mut_dict[i] = [pos,alt_at_pos]

	samples = slim.genome_dict()
	genomes = samples
	seqs = {}
	individuals_chosen = []
	while len(individuals_chosen) < number_of_sequences:
	### Add a loop here to get multiple individuals
		individual_to_choose = random.randint(1,(len(genomes.keys())/2))
		if individual_to_choose in individuals_chosen: continue
		else: individuals_chosen.append(individual_to_choose)
		seq_1 = "p1:"+str(individual_to_choose *2)
		seq_2 = "p1:"+str(individual_to_choose *2 -1)
#		brace()
		for p in [seq_1,seq_2]:
			seqs[p] = list(reference)
			for allele in genomes[p]:			
				try:
					mutation = mut_dict[int(allele)]
					seqs[p][mutation[0]] = mutation[1]
				except KeyError:pass
	return seqs,slim.name

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
parser.add_argument("--sequences","-n", 
		required = False,
		dest = "sequences", 
		type = int, 
		help = "How many individuals' sequences do you to extract?",
		default = 1)
parser.add_argument("--singletons",'-s',
                    required=False,
                    type=bool,
                    help='Do you want to include singletons? Default = False',
                    default=False)
parser.add_argument("--mut_matrix","-m", 
		required = True,
		dest = "mut_mat", 
		type = str, 
		help = "Give the name of the mutation matrix that you want to use")


args = parser.parse_args()


## Get the default output name if not supplied
if not args.output:
	if args.gz:
		output_name = args.input.replace(".txt.gz",".seq.pkl")
	else:
		output_name = args.input.replace(".txt",".seq.pkl")
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
for i in slim_iterator(args.input):
	seqs,name = generate_seq_dicts(i,mat_dict,args.sequences,singletons=args.singletons)		
	individual_dict = {}
	for j in seqs.keys():
		individual_dict[j.replace(":","_")] = "".join(seqs[j])	
	output_dict[name] = individual_dict
#pickle_jar = pickle.dump(output_dict,open(output_name,"wb"))
print "Putting all sequences into the pickle jar called", output_name



output_fasta = open(args.output,"w")
output_seq =[]

for i in output_dict.keys():
	print i
	for x in individual_dict.keys():		
		print x
		output_id = i+"_"+x+"_"
		output_fasta.write(">"+output_id+"\n"+wrap("".join(individual_dict[x]))+"\n")

