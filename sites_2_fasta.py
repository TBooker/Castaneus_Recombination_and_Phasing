#!/usr/bin/env python
import argparse,vcf, textwrap


parser = argparse.ArgumentParser(description="This script generates a fasta file, or several fasta files(LATER) from a full tabixed GZVCF and a file of sites with the haplotype info ")

parser.add_argument("combined", type = str, help="name of the combined file that has positions and haplotypes")
parser.add_argument("region", type = str, help="The chromosomal position, in tabix format, for the chromosome and position you used")
parser.add_argument("GZVCF",type = str, help = "The name of the tabixed, G-zipped VCF file that contains info on all sites") 
parser.add_argument("--output", type = str, help="The name of the output fasta, default =[haplotypes.fasta]",default = "haplotypes.fasta")

args = parser.parse_args()
chrom = args.region.split(":")[0].split("r")[1]
start = int(args.region.split(":")[1].split("-")[0])
end = int(args.region.split(":")[1].split("-")[1])

with open(args.combined) as FILE:
	for line in FILE:
		header = line
		break
full_haps = {}
hap_index = {}
for i in header.strip("\n").split("\t")[7:]:
	full_haps[i] = ""
	hap_index[header.strip("\n").split("\t").index(i)] =i 
vcf = vcf.Reader(open(args.GZVCF))

with open(args.combined) as FILE:
	index = 0
	previous_pos = start
	for line in FILE:
		if line == header:
			continue
		items = line.strip("\n").split("\t")
		var_pos = int(items[0])
		print var_pos
### The water gets a little choppy here. 
## what this section does is to grab the chunk of the genome before the 		
		if var_pos - previous_pos >1:
			chunk =vcf.fetch(chrom,start,var_pos-1)
			ref_chunk = "".join([j.REF for j in chunk])
			for key in hap_index.keys():
				full_haps[hap_index[key]]+= ref_chunk+items[key]
		elif var_pos-previous_pos ==1:
			for key in hap_index.keys():
				full_haps[hap_index[key]]+= items[key]
		previous_pos = var_pos
		index +=1

	####
### Var_pos should now represent the last line in the file...
if end - var_pos >1:
	chunk =vcf.fetch(chrom,var_pos+1,end)
	ref_chunk = "".join([j.REF for j in chunk])
	for key in hap_index.keys():
		full_haps[hap_index[key]]+= ref_chunk
elif var_pos-previous_pos ==1:
	pass
## You should now have a dictionary with n haplotypes
## You now need to write these to a FASTA file...
out_fasta = open(args.output,"w")
for i in full_haps.keys():
	out_fasta.write(">" + i + "\n" + "".join(textwrap.wrap(full_haps[i],width =70)) +"\n")



