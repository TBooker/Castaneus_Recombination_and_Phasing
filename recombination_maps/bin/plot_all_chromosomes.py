# This script generates a table that is convenient for ggplot2
# It takes each of the recombination maps and gets the scaled recombination rate at each point

import pandas as pd
import numpy as np
import sys, argparse


parser = argparse.ArgumentParser(description="Converts the recombination maps in 4Ner into cM/Mb by scaling accoridng to the Cox et al (2009) maps")

parser.add_argument("-w","--wsize", 
		required = True,
		dest = "wsize", 
		type =float, 
		help = "The window size to be used for the conversion")
parser.add_argument("-s","--slide", 
		required = True,
		dest = "slide", 
		type =float, 
		help = "The distance to slide. For plotting a smoothed line, use s=1 for looking at correlations, use s=w")

args = parser.parse_args()

cox = ("/home/booker/mouse_genome/recombination_maps/dan_files_test/Revised_HSmap_SNPs.csv")
brun = ("/home/booker/mouse_genome/recombination_maps/brunschwig_S1.csv")

x1 = pd.read_csv(cox)
x2 = pd.read_csv(brun)

first = True
diggity= pd.DataFrame()

for chrom in range(1,20)[::-1]:
	chro = str(chrom)
	print 'running chromosome:',chro
#	cast = ("~/project//6.Phasing_Castaneus+Recombination_Rates/recombination_maps/chr"+chr+"/chr"+chr+".con.map")

	casta = ("/home/booker/project/test_f_carolina/recombination_grabber/by_chrom/M.m.castneus."+chro+".map")
	print casta
	try:
		x3 = pd.read_csv(casta,header=None,sep='\t',names = ['chrom', 'type', 'method', 'pos', 'build37','dot','4Ner.bp','cum_4Ner.bp'])
	except IOError:
		print 'chromosome not found yet, dum-dum'
		continue
		
	d1 = x1.loc[x1['chr'] == int(chro)].copy()
	d1['cm_mb'] = (d1['ave_cM'] / (d1['build37']/1e6))

	d2 = x2.loc[x2['Chr'] == 'chr'+str(chro)].copy()
	d2['build37'] = (d2['Kb'] * 1e3)
	d2['ave_cM'] = np.cumsum(d2['4Ner/kb'])/sum(d2['4Ner/kb']) * max(d1['ave_cM']) # use name from Revised_HSmap_SNPs.csv

	
	d3 = x3.loc[x3['chrom'] == 'chr'+str(chro)].copy()
	d3['ave_cM'] = np.cumsum(d3['4Ner.bp'])/sum(d3['4Ner.bp']) * max(d1['ave_cM']) # use name from Revised_HSmap_SNPs.csv

	slide = args.slide
	wsize = args.wsize
	d3.to_csv(str(chro)+'.dataframe.csv')
	print d3
	continue
	s = np.arange(0, int(max(d3['build37'])/1e6), slide)
	sump = pd.DataFrame()
	sump['start'] = s
	sump['end'] = s+wsize
	sump['mid'] = (sump['start']+sump['end'])/2

	rate_dict = {}
	rate_dict['ind1'] =np.zeros(len(sump.index))
	rate_dict['ind2'] =np.zeros(len(sump.index))
	rate_dict['ind3'] =np.zeros(len(sump.index))

	for r in sump.index:
	
		ind1 = (d1['build37']/1e6 > sump['start'][r]) & (d1['build37']/1e6 <= sump['end'][r])
		try:
			rate_dict['ind1'][r] =  (d1['ave_cM'][ind1].iloc[-1] - d1['ave_cM'][ind1].iloc[0])/wsize
		except IndexError:
			continue
		ind2 = (d2['build37']/1e6 > sump['start'][r]) & (d2['build37']/1e6 <= sump['end'][r])
		try:
			rate_dict['ind2'][r] =  (d2['ave_cM'][ind2].iloc[-1] - d2['ave_cM'][ind2].iloc[0])/wsize
		except IndexError:
			continue
			
		ind3 = (d3['build37']/1e6 > sump['start'][r]) & (d3['build37']/1e6 <= sump['end'][r])
		try:
			rate_dict['ind3'][r] =  (d3['ave_cM'][ind3].iloc[-1] - d3['ave_cM'][ind3].iloc[0])/wsize
		except IndexError:
			continue
		#print rate_dict
	sump['cox'] = rate_dict['ind1'] 
	sump['brun'] = rate_dict['ind2'] 
	sump['cast'] = rate_dict['ind3'] 
	sump['chrom'] = 'Chromosome_'+str(chro)

	diggity = pd.concat([diggity, sump], axis=0)

diggity.to_csv(str(args.wsize)+'.'+str(args.slide)+'.dataframe.csv')
