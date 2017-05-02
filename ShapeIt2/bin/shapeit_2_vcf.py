## A script to convert the shapit output into VCF format which can be parsed by to Rob's script (vcf2fasta.py) Whic hI've adapted into taking phased data


import sys

if len(sys.argv) < 3:
	print "give the haplotype output (*.haps) of ShapeIt2 (assemble mode) to the program at the command line. AND give the name of the output file"
	sys.exit()

haps = open(sys.argv[1],"r")
output = open(sys.argv[2],"w")

output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tH12\tH14\tH15\tH24\tH26\tH27\tH28\tH30\tH34\tH36\n")

counter=0
for i in haps:
	counter +=1
	items = i.strip("\n").split(" ")

	try:
		c = int(items[0])
		chrom = "chr" +str(c)
	except ValueError:
		c = items[0]
		chrom = c
	pos = items[1]
	name = items[2]
	ref= items[3]
	alt =items[4]
	haplotypes = items[5:]
	genotypes = ["|".join(haplotypes[i:i+2]) for i  in range(0, len(haplotypes), 2)] 
		## This little doozy, takes every phased haplotype and crams them back into genotypes that can then be used to generate FASTA sequences, which can then, in turn, be put into LDhelmet
	#print genotypes
	outline = [chrom,pos,name,ref, alt,".","PASS",".", "GT"] + genotypes
	output.write("\t".join(outline)+"\n")
#19	3003069	3003069	G	C	10.3	PASS	.	GT	0/0	0/1	1/1	1/1	0/1	0/1	0/1	0/0	0/1	0/1
