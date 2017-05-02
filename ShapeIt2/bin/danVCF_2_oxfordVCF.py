#!/usr/bin/env python
### script to convert the SNPs-afs VCF files into ones that can be used by extractPIRs and ShapeIt2
#CHROM	POS		ID		REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097
#chr20	10000107	10000107	T	C	.	PASS	.	GT	0/0	0/0

# The above structure is what we are aiming for... Extract the information from a VCF file on lanner and use to extractPIRs
##15th September, Tom Booker

import vcf, argparse
from tom import brace
def indel(record):
	"""take a vcf record and return TRUE if the record is an indel"""
	for i in record.ALT:
		if i:
			if len(i) > 1:
				return True
		else:
			return False
	if len(record.REF) >1:
		return True
	else:
		return False



def convert_one_VCF_line(vcf_record,variants = True, bialleles = True,no_indels = True,HWE=False,minDP=0,maxDP=9999,minQUAL=0,minGQ=0):
	alt = vcf_record.ALT
	if bialleles == True:  # This is a condition that if satistifed prevents anything from being returned
		if len(alt) > 2: return
	if no_indels == True:
		if indel(vcf_record): return
	sample_genotypes = []
	if variants == True:
		gt_dict = {0:"0/0",1:"0/1",2:"1/1"}
	for i in vcf_record.samples:
		if variants ==True:
			if i.gt_type == None:continue
			sample_genotypes.append(gt_dict[i.gt_type])
		elif variants == False:
			sample_genotypes.append("0/0")
	if len(alt) ==1 or "." in alt: ## Is invariant as a reference allele is 'X' and if all sites have these then 'X' is all there is
		if variants ==True: return
		else: pass

		
		out_list = [str(vcf_record.CHROM),str(vcf_record.POS), \
		str(vcf_record.POS),str(vcf_record.REF),\
		".",str(vcf_record.QUAL),\
		"PASS",".","GT"] + sample_genotypes
		return "\t".join(out_list)

	else:	## Is variant as a variant call line must have either a non-reference allele as well as an 'X' 
		if "X" in alt:
			alt.remove("X")
		if float(vcf_record.QUAL) < minQUAL:
			return
		if HWE:
			try:
				if float(vcf_record.INFO["HWE"]) < 0.0002:return
			except KeyError:
				pass 
		null = False
		for sam in vcf_record.samples:
			if float(sam["GQ"]) < minGQ:
				null = True
			if float(sam["DP"]) < minDP or float(sam["DP"]) > maxDP:
				null = True
		if null == True:return

		alt = map(str,alt)
		out_list = [str(vcf_record.CHROM),str(vcf_record.POS), \
		str(vcf_record.POS),str(vcf_record.REF),\
		",".join(alt),str(vcf_record.QUAL),\
		"PASS",".","GT"]+sample_genotypes
		return "\t".join(out_list)


parser = argparse.ArgumentParser(description="Converts a VCF file in Dan's Format into one in the Format required by the Oxford Stat. Gen. suite of programs (Shapeit, Impute etc...)")

parser.add_argument("-v", "--vcf", \
			required = True, \
			metavar ="vcf", \
			type =str, \
			help = "The name, with full path if necessary, of the VCF file you want to convert. This program will output the VCF with the extension .ox.vcf")
parser.add_argument("-c", "--chr", \
			required = True, \
			metavar ="chr", \
			type =str, \
			help = "The chromosome that you want to extract. Please give it in the form that is present in the VCF file")
parser.add_argument("-o", "--output", \
			required = True, \
			metavar ="output", \
			type =str, \
			help = "The prefix to the output files")
parser.add_argument("--HWE", \
			required = False, \
			action="store_true", \
			help = "Do you want to reject all those sites that fail to pass a HWE test (p < 0.0002")
parser.add_argument("--minDP", \
			required = False, \
			metavar ="minDP", \
			type=int, \
			help = "what min depth (DP) do you want to take for each polymorphic site?",
			default=0)
parser.add_argument("--maxDP", \
			required = False, \
			metavar ="maxDP", \
			type=int, \
			help = "what max depth (DP) do you want to take for each polymorphic site?",
			default=9999)
parser.add_argument("--minGQ", \
			required = False, \
			metavar ="minGQ", \
			type=int, \
			help = "what min Genotype Quality (GQ) do you want to apply for each polymorphic site?",
			default=0)
parser.add_argument("--minQUAL", \
			required = False, \
			metavar ="minQUAL", \
			type=int, \
			help = "what min variant Quality (QUAL) do you want to apply for each polymorphic site?",
			default=0)






mouse_header = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H14	H15	H24	H27	H28	H30	H36	H40	H46	H62\n"


args = parser.parse_args()

newVCF = open(args.output+".ox.vcf","w")
newVCF.write(mouse_header)
reader = vcf.Reader(open(args.vcf,"r"))
start = reader.next().POS


for i in reader:
	x = convert_one_VCF_line(i,HWE=args.HWE,minDP=args.minDP,maxDP=args.maxDP,minQUAL=args.minQUAL,minGQ=args.minGQ)
	if x:
		newVCF.write(x+"\n")
newVCF.close()
