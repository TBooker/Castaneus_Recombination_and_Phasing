import pandas as pd
import argparse, sys
import numpy as np
import gzip
from itertools import izip

def find_hotspots(in_file, out, chrom, block_size, flank_size,simulation):
	if not simulation:
		d = pd.read_csv(in_file, sep=",", skiprows=0, header=None, 
			names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.500','p0.975'])
	elif simulation:
		d = pd.read_csv(in_file, sep=" ", skiprows=3,header=None ,
			names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.500','p0.975'])

	o = open(out, 'w')
	o.write('chr,block_start,block_end,flank_rate,block_rate,rate_ratio\n')

	chr_lengths = {0:500000,1:195471971,10:130694993,11:122082543,12:120129022,13:120421639,14:124902244,15:104043685,16:98207768,17:94987271,18:90702639,19:61431566,2:182113224,3:160039680,4:156508116,5:151834684,6:149736546,7:145441459,8:129401213,9:124595110}

	
	if chrom == str(0):
		chr_start = 1
	else:
		chr_start = 3000000
	if chrom == 'X':
		chr_end = 166650296
	else:
		chr_end = chr_lengths[int(chrom)]

	for block_start in range(chr_start, chr_end, int(block_size / 2.0)):
                block_end = block_start + block_size
                if block_end > chr_end:
                        block_end = chr_end
                left_flank_start = chr_start
                if (block_start - flank_size) >= left_flank_start:
                        left_flank_start = block_start - flank_size
                tmp = d[d.right_snp >= left_flank_start]
                right_flank_end = chr_end
                if (block_end + flank_size) <= right_flank_end:
                        right_flank_end = block_end + flank_size
                tmp = tmp[tmp.left_snp <= right_flank_end]
                
                chunk_rho = {'left': {'bp': 0, 'rho': 0}, 'block': {'bp': 0, 'rho': 0}, 'right': {'bp': 0, 'rho': 0}}
                chunks = {'left': [int(left_flank_start), int(block_start)], 'block': [int(block_start), int(block_end)], 'right': [int(block_end), int(right_flank_end)]}
	        for start, end, rate in izip(tmp.left_snp, tmp.right_snp, tmp.meanrho):
                        start = int(start)
                        end = int(end)

                        for chunk in chunks:
                                rate_mult = 0
                                diff = 0
                                # contained within
                                if start >= chunks[chunk][0] and end <= chunks[chunk][1]:
                                        rate_mult = (end - start) * rate
                                        diff = end - start
                                # hanging off to the left
                                elif start < chunks[chunk][0] and (end <= chunks[chunk][1] and end >= chunks[chunk][0]):
                                        rate_mult = (end - chunks[chunk][0]) * rate
                                        diff = end - chunks[chunk][0]
                           	# hanging off to the right
			   	elif (start >= chunks[chunk][0] and start <= chunks[chunk][1]) and end > chunks[chunk][1]:
                                        rate_mult = (chunks[chunk][1] - start) * rate
                                        diff = (chunks[chunk][1] - start)
				# hanging off on both sides
				elif start <= chunks[chunk][0] and end >= chunks[chunk][1]:
					diff = (chunks[chunk][1] - chunks[chunk][0])
					rate_mult = diff * rate
                                chunk_rho[chunk]['rho'] += rate_mult
                                chunk_rho[chunk]['bp'] += diff


		flank_rate = 'NA'
		block_rate = 'NA'
		diff = 'NA'
		if (chunk_rho['left']['bp'] + chunk_rho['right']['bp']) > 0:
			flank_rate = '%.5f' % ((chunk_rho['left']['rho'] + chunk_rho['right']['rho']) / float(chunk_rho['left']['bp'] + chunk_rho['right']['bp']))
		if chunk_rho['block']['bp'] > 0:
                	block_rate = '%.5f' % (chunk_rho['block']['rho'] / chunk_rho['block']['bp'])
		if block_rate != 'NA' and flank_rate != 'NA':	
			if float(flank_rate) > 0:
				diff = '%.5f' % (float(block_rate) / float(flank_rate))       

                o.write('%s,%s,%s,%s,%s,%s\n' % (chrom, block_start, block_end, flank_rate, block_rate, diff))
	o.close()
	return
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--file", required = True, help="The LDhelmet output file that you want to analyse")
	parser.add_argument("--chrom", required = False, help="the chromosome to analyse", default = 0)
	parser.add_argument("--block", required = True, help="size in bp for block size to evaluate as putative hotspot")
	parser.add_argument("--flank", required = True, help="size in bp for flank size to evaluate on either side of hotspot")	
	parser.add_argument("--output", required = True, help="The name for the output file")
	parser.add_argument("--simulation", action = "store_true", help="Use this flag if the input is simulation data (i.e. has the header and is space separated")

	args = parser.parse_args()
	chrom = args.chrom
	block_size = int(args.block)
	flank_size = int(args.flank)
	if not args.simulation and args.chrom == 0:
		print 'Give a proper chromosome, you Bozo'
		sys.exit()

	in_file = args.file
	out = args.output
	find_hotspots(in_file, out, chrom, block_size, flank_size, args.simulation) 

if __name__ == "__main__":
    main()
