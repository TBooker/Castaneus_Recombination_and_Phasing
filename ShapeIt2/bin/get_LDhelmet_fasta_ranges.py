import argparse, sys,gzip
### This script gets files needed for LDhelmet. The prblem with LDhelmet is that you need the computation requires a HUGE amount of memory.
### Running this for a whole chromosome is untenable, so break into chunks. The authors of LDhelmet did the very same thing in their paper




parser = argparse.ArgumentParser(description="Take a VCF file of phased haplotypes and produces a list of tabix regions for other scripts to use", usage="get_LDhelmet_fasta.py --vcf [-v] phased_sequences.vcf --chrom [-c] CHROM")

parser.add_argument("-v", "--vcf",
        required = True,
	type = str,
	help = "This is VCF file that you want to get coordinates from, NOT BGZIPPED!!!")

parser.add_argument("-c", "--chrom",
        required = True,
	type = str,
	help = "The chromosome you are using")
parser.add_argument("-o", "--output",
        required = False,
	type = str,
	help = "The name of the file that contain the TABIX format ranges",
	default = "tabix_ranges.txt")
parser.add_argument("--snps",
        required = False,
	type = int,
	help = "The number of SNPs to a fasta",
	default = 4000)
parser.add_argument("--overhang",
        required = False,
	type = int,
	help = "The number of SNPs to overhang",
	default = 200)

args = parser.parse_args()

vcf = gzip.open(args.vcf,"r")

length = 0
coords = [] ## Iterate through VCF and produce a list of coordinates
new=False
for i in vcf:
	if i.startswith("#"):continue
	length +=1
print "there are "+str(length)+" SNPs in " + args.vcf
vcf.close()
width = args.snps
LD_decay = args.overhang
coord_range =  range(1,length,width)
ranges = [[1,width+(LD_decay)]]+[[i-LD_decay,i+width+LD_decay] for i in coord_range[1:]]
print "giving a total of " + str(len(ranges)) + " ranges to get fasta files for"
start_index = 0
end_index =0
counter =0
starts =[]
ends =[]
#print ranges
######
#sys.exit()
vcf = gzip.open(args.vcf,"r")
if ranges[-1][1] > length:
	ranges[-1][1] = length
for i in vcf:
#	if range_index >len(ranges)-1:continue 
	current_start =ranges[start_index][0]
	current_end =ranges[end_index][1]
	if i.startswith("#"):continue
	counter +=1
	x = i.split("\t")
	if args.chrom != x[0]:
		print "The chromosome you speficied does not match the VCF file you gave"
		print "you gave: "+args.chrom
		print "The VCF lists: "+x[0]
		break
	pos = x[1]
	if counter == current_start:
		starts.append(pos)
		start_index += 1
		if start_index == len(ranges):	
			start_index -=1
	elif counter == current_end:
		ends.append(pos)
		end_index += 1
		if end_index == len(ranges):	
			end_index +=1
vcf.close()

ranges_out = open(args.output,"w")
for i,j,k in zip(ranges,starts,ends):

	ranges_out.write(args.chrom +":"+ j + "-" + k + "\n")

ranges_out.close()
