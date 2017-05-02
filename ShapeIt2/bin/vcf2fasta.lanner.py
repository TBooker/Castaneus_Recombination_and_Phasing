#! /usr/bin/python

import vcf, sys,re,argparse
from Bio import SeqIO
from collections import OrderedDict
#from optimizer_functions import  hap2dip_converter
def wrap(to_wrap, width=70):
    return '\n'.join(to_wrap[i:i+width] for i in range(0, len(to_wrap), width))


def hap2dip_converter(bases):
	converter = {'AA': 'A', 'TT':'T', 'GG':'G', 'CC':'C', 'AT':'W', 'GC':'S', 'AG':'R', 'CT':'Y', 'GT': 'K', 'AC': 'M', 'AN':'A','CN':'C','GN':'G','TN':'T','TA':'W', 'CG':'S', 'GA':'R', 'TC':'Y', 'TG': 'K', 'CA': 'M', 'NA':'A','NC':'C','NG':'G','NT':'T', 'NN':'N'}
	hap_base = converter[bases]
	return hap_base

# pseudo code
# for each position
# 	decide what the alleles are
# 	if its invariant
# 		apply ref allele
# 	if its variant
# 		for each sample
# 			decide its genotype
# 			insert its alleles into the sequence
# 				if insertion change position to position + insertion
# 				if deletion remove the next N bases
# 				if SNP 
# 		if consensus
# 			decide consensus genotype
# 			insert into consensus
# def consensus_position(vcf_record, position,purity=0.5, min_GQ=0, indels=False, samples='all'):
# 	#################################################################
# 	""" takes a vcf_file and position and returns the majority allele """
# 	#################################################################
# 	if vcf_record.POS != position:
# 		consensus_allele = 'N'
# 		return consensus_allele
# 	if indel(vcf_record) or vcf_record.is_snp:
# 		#there are some decisions to make
# 		GTs = [i['GT'] for i in vcf_record.samples if i['GQ']>min_GQ ]
# 		alleles = re.findall('[0-9]',"".join(GTs)) #so this is a bunch of numbers 
# 		majority = Counter(re.findall('[0-9]',"".join(GTs))).most_common(1)[0] #this will be tuple (allele[str], count[int]) most common number
# 		#print vcf_record.POS, Counter(re.findall('[0-9]',"".join(GTs)))
# 		if majority[1]/float(len(alleles)) > purity: #if its sufficiently common
# 			maj_allele = majority[0]
# 			if maj_allele == '0':
# 				consensus_allele = str(vcf_record.REF)
# 				if indels == False:
# 					consensus_allele = str(vcf_record.REF)[0]
# 			else:
# 				consensus_allele=str(vcf_record.ALT[int(maj_allele)-1])
# 				if indels == False and len(consensus_allele) > 1:
# 					consensus_allele = str(vcf_record.REF)[0]
# 		else:
# 			consensus_allele=str(vcf_record.REF)
# 			if indels == False:
# 				consensus_allele = str(vcf_record.REF)[0]
# 	else:
# 		consensus_allele=str(vcf_record.REF)
# 		if indels == False:
# 			consensus_allele = str(vcf_record.REF)[0]
# 	return consensus_allele 
# 	#it must be the minor allele unless there are two genotypes (ie 50:50) or many alt alleles
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
def rc(sequence):
	complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'R':'Y', 'Y':'R', 'K':'M', 'M':'K' ,'W':'W', 'S':'S', '?': '?', 'X':'X', '.': '.', '-':'-', 'B':'V', 'D':'H', 'H':'D', 'V':'B'}
	return ''.join([complements[base] for base in sequence][::-1])
def get_genotype_calls(record, samples,min_GQ=0,min_DP=0,min_QUAL=0,min_MQ=0, consensus=False, sample_ploidy=2, vcf_ploidy=2, fill_w_reference=False, missing = 'N', filter_heterozygote_haploids =True, parse_indels = False,phased =False):
	invariant = False
	variant = False
	chromosome = record.CHROM
	#now if it were a new gene everything is empty otherwise it just carries on
	position = record.POS
	#make allele dict
	if not record.ALT[0]:
		allele_dict = {0:record.REF}
		invariant = True
		#invariant
	elif record.ALT[0]: #variant
		variant = True
		allele_dict = {0:record.REF}
		idx=1
		for alt in record.ALT: 
			allele_dict[idx] = str(record.ALT[idx-1])
			idx+=1
	#now make genotype calls
	site_genotypes = OrderedDict()
	if record.QUAL == None:
		record.QUAL = 1
	if record.QUAL > min_QUAL:# and record.INFO['MQ']> min_MQ: #minimum quality site
		# if consensus == True:
		# 	consensus_allele = consensus_position(record, position, min_GQ=min_GQ, samples=samples)
		for s in samples:
			if invariant:
				#gt, dp = record.genotype(s)['GT'], record.genotype(s)['DP']
				gt = record.genotype(s)['GT']
				if vcf_ploidy == 2:
					if sample_ploidy==2 and phased == False:
						if gt and dp>min_DP:
							site_genotypes[s+"_1"] = (allele_dict[int(gt.split("|")[0])])
							site_genotypes[s+"_2"] = (allele_dict[int(gt.split("|")[1])])
						else:#uncallable
							if fill_w_reference: 
								site_genotypes[s+"_1"] = (record.REF, record.REF) 
								site_genotypes[s+"_2"] = (record.REF, record.REF) 

							else: 
								site_genotypes[s+"_1"] = (missing)
								site_genotypes[s+"_2"] = (missing, missing)
					elif sample_ploidy ==1 or phased == True:
						if gt and phased == False:
							#print dp
							# filter heterozygotes?
							if gt.split("/")[0] != gt.split("/")[1]:
								if filter_heterozygote_haploids: #heterozygote
									site_genotypes[s+"_1"] = (missing)
									site_genotypes[s+"_2"] = (missing)

								else:
									sys.stderr.write("I don't have a good way of dealing with haploids that appear heterozygous")
									sys.exit()
							else: #homozygote appearance
								site_genotypes[s] = (allele_dict[int(gt.split("/")[0])],)
						if gt and phased == True:
							if gt.split("|")[0] != gt.split("|")[1]:
									site_genotypes[s+"_1"] = (record.REF)
									site_genotypes[s+"_2"] = (record.REF)
							else: #homozygote appearance
								site_genotypes[s] = (allele_dict[int(gt.split("|")[0])],)

						else: #no call
							if fill_w_reference: site_genotypes[s] = (record.REF,) 
							else: site_genotypes[s] = (missing,)
					else:
						sys.stderr.write("we don't deal with your weird ploidy")
						sys.exit()
				elif vcf_ploidy==1:
					if gt and dp>min_DP:
						 site_genotypes[s] = (allele_dict[int(gt.split("/")[0])],)
					else:#uncallable
						if fill_w_reference: site_genotypes[s] = (record.REF,) 
						else: site_genotypes[s] = (missing,)
			if variant:
				
				#gt, dp, gq = record.genotype(s)['GT'], record.genotype(s)['DP'], record.genotype(s)['GQ']
				#print record
				#print s
				gt = record.genotype(s)['GT']#, record.genotype(s)['DP'], record.genotype(s)['GQ']

				if vcf_ploidy == 2:
					if sample_ploidy==2 and phased == False:
						if gt:

#						if gt and dp>min_DP and gq>min_GQ:
							site_genotypes[s+"_1"] = (allele_dict[int(gt.split("|")[0])])
							site_genotypes[s+"_2"] = (allele_dict[int(gt.split("|")[1])])
							

						else:#uncallable
							if fill_w_reference: site_genotypes[s] = (record.REF, record.REF) 
							else: site_genotypes[s] = (missing, missing)
					elif sample_ploidy ==1 or phased == True:
						if gt:
#						if gt and dp>min_DP and gq>min_GQ:
							# filter heterozygotes?
							if phased == True:
								site_genotypes[s+"_1"] = (allele_dict[int(gt.split("|")[0])]).lower()
								site_genotypes[s+"_2"] = (allele_dict[int(gt.split("|")[1])]).lower()
							if sample_ploidy ==1:
								if filter_heterozygote_haploids and gt.split("/")[0] != gt.split("/")[1]:
									site_genotypes[s] = (missing,)
								elif not filter_heterozygote_haploids and gt.split("/")[0] != gt.split("/")[1]:
									sys.stderr.write("I don't have a good way of dealing with haploids that appear heterozygous")
									sys.exit()
								else: #homozygote appearance
									site_genotypes[s] = (allele_dict[int(gt.split("|")[0])],)
						else: #no call
							if fill_w_reference:
								if phased == False: 
									site_genotypes[s] = (record.REF,) 
								if phased == True: 
									site_genotypes[s+"_1"] = (record.REF)
									site_genotypes[s+"_2"] = (record.REF)  
	
							else:
								site_genotypes[s] = (missing,)
					else:
						sys.stderr.write("we don't deal with your weird ploidy")
						sys.exit()
				elif vcf_ploidy==1:
					if gt and dp>min_DP:
						 site_genotypes[s] = (allele_dict[int(gt.split("/")[0])],)
					else:#uncallable
						if fill_w_reference: site_genotypes[s] = (record.REF,) 


						else: site_genotypes[s] = (missing,)
	else: #this site is low qual - what do we do?
		for s in samples:
			if sample_ploidy==2:
				if fill_w_reference: site_genotypes[s] = (record.REF, record.REF) 
				else: site_genotypes[s] = (missing, missing)
			elif sample_ploidy==1:
				if fill_w_reference: site_genotypes[s] = (record.REF,) 
				else: site_genotypes[s] = (missing,)
	return site_genotypes


def vcf2fasta(vcf_file, chromosome, start, end, reference_fasta, min_GQ=0,min_DP=0,min_QUAL=0,min_MQ=0,samples='all', consensus=False, sample_ploidy=2, vcf_ploidy=2, fill_w_reference=False, missing = 'N', filter_heterozygote_haploids =True, filter_heterozygous_haploid_sites =True, parse_indels = False, phased = False):
	#things to think about
	# if parse_indels and sample_ploidy ==1 and vcf_ploidy==1 and filter_heterozygote_haploids ==False <- this is an impossible set of situations 
	vr = vcf.Reader(filename=vcf_file)
	current_chromosome = chromosome
	skipable = []
	first_position=True 
	ref_dict = SeqIO.to_dict(SeqIO.parse(open(reference_fasta, 'r'), 'fasta'))
	output_seqs = OrderedDict()
	if samples[0] == 'all': samples = vr.samples
	for s in samples:
		if phased==True:
			output_seqs[s+"_1"] = list(ref_dict[chromosome].seq[start-1:end])
			output_seqs[s+"_2"] = list(ref_dict[chromosome].seq[start-1:end])
		else:
			output_seqs[s] = list(ref_dict[chromosome].seq[start-1:end])
	#output_seqs['reference'] = str(ref_dict[chromosome].seq[start-1:end])
	if consensus: consensus_seq = list(ref_dict[chromosome].seq[start-1:end])
	previous_position = start-1
	for record in vr.fetch(chromosome, start, end):

		position = record.POS ## this is for a fuck up in the getting of chromosome names
		if position == previous_position +1 and not record.ALT[0] and fill_w_reference:
			pass #if its invariant and fill_w_reference - there is nothing to do so skip the rest
		else:
			site_genotypes = get_genotype_calls(record, samples,min_GQ=min_GQ,min_DP=min_DP,min_QUAL=min_QUAL,min_MQ=min_MQ, consensus=consensus, sample_ploidy=sample_ploidy, vcf_ploidy=vcf_ploidy, fill_w_reference=fill_w_reference, missing = missing, filter_heterozygote_haploids =filter_heterozygote_haploids, parse_indels = parse_indels,phased = phased)
			if first_position and position != start: #fill in the 5' end
				if not fill_w_reference:
					#swap in missing symbol for all the missing bases
					for s in output_seqs:
						for base in range(0,position-start):
							output_seqs[s][base] = missing
				else: pass #this is assuming the reference is all parsed into output seqs
			first_position = False
			#### Missing reference data #################################################################################################
			if position > previous_position + 1:
				
				#fill in with N's, ref or if it was deleted leave it alone
				for s in output_seqs.keys():
					if fill_w_reference: 
						pass
					else:
						for i in range((previous_position-start)+1, (previous_position-start)+(position-previous_position)):

							output_seqs[s][position-start] = missing
				# now if we set previous position to position-1 it will carry on normally
				previous_position = position-1
			#### Normal Position #################################################################################################
			if position == previous_position + 1: #normal_position
				if record.ALT and fill_w_reference and phased == False: #
					pass
				else:
					for s in output_seqs.keys():
						
						#if the site is deleted ("-") then skip it
						if sample_ploidy == 1 and phased == False:
							if filter_heterozygous_haploid_sites and 	record.num_het > 0:
								output_seqs[s][position-start] = missing
							else:
								output_seqs[s][position-start] = site_genotypes[s]
						elif sample_ploidy == 2 and phased ==True:
							output_seqs[s][position-start] = site_genotypes[s]			## This makes haploid sequences for a phased diploid
						elif sample_ploidy == 2 and phased ==False:
							#this currently is for making just one diploid sequence ie with ambiguity calls
							output_seqs[s][position-start] = hap2dip_converter("".join(site_genotypes[s]))
							# pseudocode
							# if sample == homozygous:
							# 	output_seqs[s][position-start] = site_genotypes[s][0]
							# elif sample == heterozygous:
							# 	if using 2 seqs:
							# 		output_seqs[s][0][position-start] = site_genotypes[s][0]
							# 		output_seqs[s][1][position-start] = site_genotypes[s][1]
							# 	elif using 1 diploid seq:
							# 		output_seqs[s][position-start] = hap2dip(site_genotypes[s])
			####   INDEL   ##################################################################################################################
			elif position == previous_position:
				#print "indel?", site_genotypes
				if parse_indels == True:
					#we are in the dreaded INDEL!!!
					if not record.ALT[0]: #invariant????
						pass #forget it and move on!!!
					elif sample_ploidy == 1 or phased == True:
						#if there are heterozygotes skip site
						#now there should only be homozygotes
						### insertion ###
						if len(record.REF) < max(len(i) for i in record.ALT):
							if len(record.REF) > 1:
								sys.stderr.write("Not equipped to deal with these complex indels" + \
												 "\t".join([str(x) for x in [record.CHROM, record.POS, record.REF, ",".join([str(i) for i in record.ALT])]]) + "\n") 
							else:
								#alter the base at position to match the insertion
								#get len of alleles, pad the non-inserted individuals with a string of ----'s after the allele
								alleles = [i[0] for i in set(site_genotypes.values())]
								idict={}
								for a in alleles:
									idict[a] = str(a) + "-"*(max(len(i) for i in alleles)-len(a))
								for s in samples:
									output_seqs[s][position-start] = idict[site_genotypes[s][0]]
						### deletion ###
						if len(record.REF) > max(len(i) for i in record.ALT):
							if max(len(i) for i in record.ALT) > 1:
								sys.stderr.write("Not equipped to deal with these complex indels" + "\t".join([str(x) for x in [record.CHROM, record.POS, record.REF, ",".join([str(i) for i in record.ALT])]]) + "\n") 
							else:
								# the number of deleted bases needs to be removed from the following positions
								# The current position is left and len(ALT) - len(REF) deletions are made to the following bases
								for s in samples:
									del_length = len(record.REF) - len(site_genotypes[s][0])
									for i in range (1,1+del_length):
										idx = (position-start) + i
										output_seqs[s][idx] = "-"
					elif sample_ploidy ==2:
						sys.stderr.write("I haven't written the code for this")
						sys.exit()
						# pseudocode
						# if its a heterozygote indel we need two sequences to represent it.
						# otherwise a flag needs to go up maybe lower case for hemizygous bases.
				else:pass
			
			####   something is weird ie position isn't indel isn't missing and isn't normal...VCF is broken? ################################
			else:
				pass
				#print "ORPHAN:"
		previous_position = position
	
	return output_seqs



parser = argparse.ArgumentParser(description="Export FASTA format sequences from a VCF", usage="vcf2fasta.py [options] my_vcf.vcf")
parser.add_argument("-v", "--vcf_file",
					required = True,
					type=str,
					help="This is the vcf input file. It must tabix indexed [required]")
parser.add_argument("-r", "--reference",
					required = True,
					type=str,
					help="This is the genome reference file used to make the vcf. Only the chromosome overlapping the region of interest is required [required]")
parser.add_argument("-i", "--regions",
					required = True,
					type=str,
					nargs="+",
					help="This is the regions of the vcf to generate a FASTA for. In SamTools format chromosome:start-end [required]")
parser.add_argument("-s", "--samples",
					required = False,
					default='all',
					type=str,
					nargs="+",
					help="This is the region of the vcf to generate a FASTA for. In SamTools format chromosome:start-end [required]")
parser.add_argument("--consensus",
					action="store_true",
					required=False,
					help="If invoked a consensus will be generated [False]")
parser.add_argument("--concatenate",
					action="store_true",
					required=False,
					help="If invoked this will patch all the sequences derived from the regions end 2 end - useful for CDS [False]")
parser.add_argument("--reverse_complement",
					action="store_true",
					required=False,
					help="If invoked this will reverse complement the extracted sequence useful for CDS on - strand [False]")
parser.add_argument("--sample_ploidy",
					default=1,
					type=int,
					help="The ploidy of the sample [1]")
parser.add_argument("--vcf_ploidy",
					default=2,
					type=int,
					help="The ploidy assumed when the vcf was called [2]")
parser.add_argument("--fill_w_reference",
					action="store_true",
					required=False,
					help="If invoked missing and low quality sites will be assumed to match the reference [False]")
parser.add_argument("--filter_heterozygote_haploids",
					action="store_false",
					required=False,
					help="If True apparently heterozygous genotypes in haploids will be assumed to be errors and marked as missing or low quality [True]")
parser.add_argument("--filter_heterozygote_haploid_sites",
					action="store_false",
					required=False,
					help="If True if any individual is apparently heterozygous the whole site in all individuals is assumed to be an error and marked as missing or low quality [True]")
parser.add_argument("--missing",
					required = False,
					type=str,
					default = "N",
					help="This character will be used to denote missing or low quality sites [N]")
parser.add_argument("-g", "--min_GQ",
					dest="min_GQ",
					default=0,
					type=int,
					help="This is the minimum GQ score that a variant site can have [0]")
parser.add_argument("--min_QUAL",
					dest="min_QUAL",
					default=0,
					type=int,
					help="This is the minimum QUAL score that a site can have to be considered [0]")
parser.add_argument("--min_MQ",
					dest="min_MQ",
					default=0,
					type=int,
					help="This is the minimum MQ score that a site can have to be considered [0]")
parser.add_argument("--min_DP",
					dest="min_DP",
					default=0,
					type=int,
					help="This is the minimum depth (DP) that a site can have to be considered [0]")
parser.add_argument("--parse_indels",
					action="store_true",
					required=False,
					help="If invoked indels will be called and inserted into sequences. Note that this will mean the annotations of genes may be altered [False]")
parser.add_argument("--phase_bool",
					action="store_true",
					required=False,
					help="If invoked sample genotpyes will be assumed to be separated by a pipe '|' and you will get the phased haplotypes out")

#############
args = parser.parse_args()
vcf_file = args.vcf_file
reference = args.reference
regions = args.regions
samples = args.samples
consensus = args.consensus
concatenate = args.concatenate
reverse_complement = args.reverse_complement
sample_ploidy = args.sample_ploidy
vcf_ploidy = args.vcf_ploidy
fill_w_reference = args.fill_w_reference
filter_heterozygote_haploids = args.filter_heterozygote_haploids
filter_heterozygote_haploid_sites = args.filter_heterozygote_haploid_sites
missing = args.missing
min_GQ = args.min_GQ
min_QUAL = args.min_QUAL
min_MQ = args.min_MQ
min_DP = args.min_DP
parse_indels = args.parse_indels
phase_bool = args.phase_bool



# print 'vcf_file', vcf_file
# print 'reference', reference
# print 'regions', regions
# print 'samples', samples
# print 'consensus', consensus
# print 'concatenate', concatenate
# print 'reverse_complement', reverse_complement
# print 'sample_ploidy', sample_ploidy
# print 'vcf_ploidy', vcf_ploidy
# print 'fill_w_reference', fill_w_reference
# print 'filter_heterozygote_haploids', filter_heterozygote_haploids
# print 'filter_heterozygote_haploid_sites', filter_heterozygote_haploid_sites
# print 'missing', missing
# print 'min_GQ', min_GQ
# print 'min_QUAL', min_QUAL
# print 'min_MQ', min_MQ
# print 'min_DP', min_DP
# print 'parse_indels', parse_indels


# common samples
# samples = 'all'
# samples = ['CC3060',  'CC3064',  'CC3065',  'CC3068',  'CC3069',  'CC3071',  'CC3072',  'CC3076',  'CC3078',  'CC3086'] #quebec mt+
# samples = ['CC3059', 'CC3061', 'CC3062', 'CC3063', 'CC3073', 'CC3075', 'CC3079', 'CC3082', 'CC3083', 'CC3084'] #quebec mt-
# samples = ['CC124', 'CC1373', 'CC1690', 'CC1952', 'CC2342', 'CC2343', 'CC2344', 'CC2931', 'CC2932', 'CC2935', 'CC2936', 'CC2937', 'CC2938'] #all wt species wide
# samples = ['CC124','CC1952', 'CC2938','CC2342', 'CC2935','CC2931' ] #wild mt-
# samples = ['CC2343','CC2936', 'CC1373','CC2344', 'CC2937','CC1690', 'CC2932'] #wild mt+
# reference_fasta = '/home2/data/genomes/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta'

output_seqs = OrderedDict()
for region in regions:
	chromosome, start, end = re.split("[-:]",region)[0], int(re.split("[-:]",region)[1]), int(re.split("[-:]",region)[2])
	seqs = vcf2fasta(vcf_file, chromosome, start, end, reference, \
	min_GQ=min_GQ,min_DP=min_DP,min_QUAL=min_QUAL,min_MQ=min_MQ,samples=samples, consensus=consensus, \
	sample_ploidy=sample_ploidy, vcf_ploidy=vcf_ploidy, fill_w_reference=fill_w_reference, missing = missing, \
	filter_heterozygote_haploids =filter_heterozygote_haploids, filter_heterozygous_haploid_sites = filter_heterozygote_haploid_sites, parse_indels = parse_indels,phased = phase_bool)
	if concatenate:
		for s in samples:
			if s not in output_seqs: output_seqs[s]=seqs[s]
			else: output_seqs[s] += seqs[s]
	else:
		for s in seqs:
			if reverse_complement: sys.stdout.write(">%s|%s\n%s\n" %(s,region,wrap(rc("".join(seqs[s]).upper()))))
			else: sys.stdout.write(">%s|%s\n%s\n" %(s,region,wrap("".join(seqs[s]).upper())))

if concatenate: #### REGION IS WRONG
	coords = []
	for region in regions:
		chromosome, start, end = re.split("[-:]",region)[0], int(re.split("[-:]",region)[1]), int(re.split("[-:]",region)[2]) 
		coords+=[start,end]
	region="%s:%i-%i" %(chromosome,min(coords),max(coords))
	for s in output_seqs:
		if reverse_complement: sys.stdout.write(">%s|%s\n%s\n" %(s,region,wrap(rc("".join(output_seqs[s]).upper()))))
		else: sys.stdout.write(">%s|%s\n%s\n" %(s,region,wrap("".join(output_seqs[s]).upper())))



"""
eg deletion (quebec_wt.vcf.gz)
19083 A [None]
19083 ACATCTACGAGCC [A]
1/1, None, None, None, None, None, None, None, None, 1/1, None, 0/1, 1/1, None, 0/0, 0/0, 0/0, 0/0, 0/0, 0/0, 


eg insertion (quebec_wt.vcf.gz)
get_genotype_calls(insertion, samples,min_GQ=0,min_DP=0,min_QUAL=0,min_MQ=0, consensus=False, sample_ploidy=1, vcf_ploidy=2, fill_w_reference=False, missing = 'N', filter_heterozygote_haploids =True, parse_indels = True)

19097 A [None]
19097 A [AGGGTATGTCGGAACGCAAGCGATC, AAAAGGCAA]
1/2, None, 2/2, 2/2, 2/2, 2/2, 2/2, 2/2, 2/2, 1/1, None, 2/2, 2/2, 2/2, 0/0, 2/2, 0/0, 0/0, 0/0, 2/2, 

>>> alleles = set(site_genotypes.values())
>>> idict={}
>>> for a in alleles:
...     idict[a] = str(a) + "-"*(max(len(i) for i in alleles)-len(a))
... 
>>> 
>>> for i in idict:
...     print idict[i], i
... 
A------------------------ A
AAAAGGCAA---------------- AAAAGGCAA
N------------------------ N
AGGGTATGTCGGAACGCAAGCGATC AGGGTATGTCGGAACGCAAGCGATC

# 1234567
# 0123456
# 4567890
# ABCDEFG
# start = 4
# pp =6
# position =9
# missing
# from pp 6 
# 7,8 are missing
# pp-start  to (pp-start)+(pos-pp)
# 6-4 to (6-4)+(9-6)
# 2:5
# deletion
position = 6 
7 and 8 are deleted
so REF = 678 ALT = 6
for i in range (1,1+del_length)
	idx = (position-start) + i
"""

