import sys,gzip

if len(sys.argv) < 3:
	print "Give the name of the .haps.gz file and the .ox.vcf.gz file at the command line\npython\tX.haps.gz\tX.ox.vcf"
	sys.exit()

haplo_file = gzip.open(sys.argv[1])
vcf_file = gzip.open(sys.argv[2])
count = 0

vcf_dict = {}
hap_dict = {}

ind = "H40"
reverse = False
het_count = 0
error_count = 0
errors = []
hap_1 =""
hap_2 =""
true_1=""
true_2=""
for h,v in zip(haplo_file,vcf_file):	
	if count==0:
		print h
		for h_key in h.strip("\n").strip("#").split(" "):
			hap_dict[h_key] = h.strip("\n").strip("#").split(" ").index(h_key)
		for v_key in v.strip("\n").strip("#").split("\t"):
			#print v_key
			vcf_dict[v_key] = v.strip("\n").strip("#").split("\t").index(v_key)
		count +=1
		continue
## Can filter out variants from the VCF

	vcf = v.strip("\n").strip("#").split("\t")

	## For the VCF called by Dan, take only homozygous variants as the truth 

	count +=1

	hap = h.strip("\n").strip("#").split(" ")


	alleles = {0:hap[hap_dict["REF"]],1:hap[hap_dict["ALT"]]}
	true_base_1 = int(vcf[vcf_dict[ind]].split("/")[0])
	true_base_2 = int(vcf[vcf_dict[ind]].split("/")[1])

	if reverse == False:
		inf_base_1 = int(hap[hap_dict[ind+"_1"]])
		inf_base_2 = int(hap[hap_dict[ind+"_2"]])
	if reverse == True:
		inf_base_1 = int(hap[hap_dict[ind+"_2"]])
		inf_base_2 = int(hap[hap_dict[ind+"_1"]])
	
	if true_base_1 != true_base_2:
		if true_base_1 != inf_base_1:
		#	print "SWITCH"
			errors.append(count)
			error_count +=1
			het_count += 1
			reverse = True
		else:
			het_count +=1
	hap_1 += str(inf_base_1)
	hap_2 += str(inf_base_2)
	true_1 += str(true_base_1)
	true_2 += str(true_base_2)
	if count == 150:
 		break

for o,i,j,k,l in zip(xrange(len(hap_1)),hap_1,hap_2,true_1,true_2):
	if o+2 not in errors:
		print o," ",i,j," ",k,l
	else: 
		print o," ",i,j," ",k,l, "< switch"

print float(error_count)/het_count

print het_count,error_count

