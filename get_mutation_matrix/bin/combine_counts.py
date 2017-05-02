import glob

big_dict = {}

for i in glob.glob("*_mut_rates.txt"):
	z = open(i,"r").readlines()
	if len(big_dict.keys()) == 0:
		for x in z:
			y = x.strip("\n").split("\t") 
			mut = ":".join([y[0],y[1]])
			big_dict[mut] = int(y[2])
	else:
		for x in z:
			y = x.strip("\n").split("\t") 
			mut = ":".join([y[0],y[1]])
			big_dict[mut] += int(y[2])
for i in big_dict.keys():
	print i,big_dict[i]
M_list = []
for i in "ACGT":
	f_i = big_dict[i+":."]
	M_temp = []
	for j in "ACGT":
		if j == i:continue
		M_temp.append(float(big_dict[i+":"+j]))
	M_list.append(sum(M_temp)/f_i)
M = max(M_list)

matrix = {}
for i in "ACGT":
	f_i = big_dict[i+":."]
	matrix[i] = {}
	for j in "ACGT":
		if j != i:
			f_ij = big_dict[i+":"+j]
			matrix[i][j] = float(f_ij)/f_i/M
		else: continue


for i in "ACGT":
	matrix[i][i] = 1- sum([matrix[i][k] for  k in "ACGT" if k!=i])

output = open("rate_matrix.txt","w")
output.write("\t".join([".","A","C","G","T"])+"\n")
for i in "ACGT":
	line = []
	for l in "ACGT":
		line.append(str(round(matrix[i][l],2)))
	output.write(i+"\t"+"\t".join(line)+"\n")
output.close()

