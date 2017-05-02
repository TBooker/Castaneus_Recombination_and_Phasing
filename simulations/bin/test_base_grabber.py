import sys, random

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
#	while alt == original:
#		alt = random.choice("ATCG")
#	return alt


test_dict = {}
for o in "ACGT":
	test_dict[o] = 0
mut_mat = {}

for k in open(sys.argv[1],"r"):
	cc = k.strip("\n").split(" ")
	mut_mat[cc[0]] = float(cc[1])

#lengths = [mut_mat[i+":."] for i in "ACGT"]
#proportions = [k/sum(lengths) for k in lengths]
#for kk in xrange(1000000): 	
#	random_position = random.uniform(0,1)
#	accumulator = 0
#	for i,j in zip(proportions,"ACGT"):
#		accumulator += i
#		if random_position > accumulator - i and random_position < accumulator :
#			base_selected = j
#	test_dict[base_selected]+=1

#for i in test_dict.keys():
#	print i,test_dict[i]

mut_mat_keys = [i for i in mut_mat.keys() if i[2] != "."]
mut_freqs = {}
mut_counts = {}
for l in "ACGT":
	mut_counts[l] = sum([mut_mat[l+":"+p] for p in  "ACGT"])
for l in "ACGT":
	for k in "ACGT":
		mut_freqs[l+":"+k] = mut_mat[l+":"+k]/mut_counts[l]


#for i in mut_mat:
#	print i,mut_mat[i]
#for i in mut_counts:
#	print i,mut_counts[i]
#for i in mut_freqs:
#	print i,mut_freqs[i]
for i in mut_freqs:
	print i,mut_freqs[i]
bases = []
for i in xrange(1000):
	bases.append(mutate("A",mut_freqs))

print "A",bases.count("A")	
print "C",bases.count("C")	
print "G",bases.count("G")	
print "T",bases.count("T")	


