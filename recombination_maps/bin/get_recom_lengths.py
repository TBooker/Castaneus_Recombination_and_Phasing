import sys
chromosome_lengths = {1:197195432,10:130694993,11:122082543,12:120129022,13:120421639,14:124902244,15:104043685,16:98207768,17:94987271,18:90702639,19:61431566,2:182113224,3:160039680,4:156508116,5:151834684,6:149736546,7:145441459,8:129401213,9:124595110}
chromo_rates = {1:0.61e-8, 2:0.54e-8, 3:0.44e-8, 4:0.54e-8, 5:0.56e-8, 6:0.45e-8, 7:0.52e-8, 8:0.61e-8, 9:0.58e-8, 10:0.61e-8, 11:0.68e-8 ,12:0.51e-8, 13:0.52e-8, 14:0.64e-8, 15:0.63e-8, 16:0.54e-8, 17:0.57e-8, 18:0.46e-8, 19:1.04e-8}
output = open("chromosome_length_summary.csv","w")

for c in range(1,20): 
	previous_pos = 3000000
	chro = "chr"+str(c)
	rho_sum = 0 
	bases = 0
	for i in open(sys.argv[1]):
		x = i.split("\t")
		if x[0] == "chrom":continue
		if x[0] == chro:
#			print x
#			print previous_pos
			
			rho = float(x[6])
			pos = int(x[3])
	
			bases+= (pos-previous_pos)
			rho_sum +=(pos-previous_pos)*rho		
			previous_pos = pos	
	bases+=(chromosome_lengths[c]-previous_pos)
	rho_sum +=(chromosome_lengths[c]-previous_pos)*rho
	d_rho_sum = 0 
	d_bases = 0
	d_previous_pos = 3000000

	outlist_c = map(str,["castaneus",c,bases,rho_sum,((rho_sum/bases)/chromo_rates[c])/4])
	output.write(",".join(outlist_c)+"\n")


	for i in open("/home/booker/mouse_genome/recombination_maps/brunschwig_S1.csv"):
		x = i.strip("\r\n").split(",")
		if x[0]!=chro:continue
		else:
			d_pos = float(x[1])*1e3
			d_rho = float(x[2])/1e3
			d_bases+= (d_pos-d_previous_pos)
			d_rho_sum +=(d_pos-d_previous_pos)*d_rho		
			d_previous_pos = d_pos
	d_bases+=(chromosome_lengths[c]-d_previous_pos)
	d_rho_sum +=(chromosome_lengths[c]-d_previous_pos)*d_rho

	#print
	#print "castaneus"
	#print bases
	#print rho_sum	
	#print ((rho_sum/bases)/6.1e-9)/4

	print
	print "domesticus",c
	print d_bases
	print d_rho_sum	
	print ((d_rho_sum/d_bases)/chromo_rates[c])/4
	outlist_d = map(str,["brunschwig",c,d_bases,d_rho_sum,((d_rho_sum/d_bases)/chromo_rates[c])/4])
	output.write(",".join(outlist_d)+"\n")
#	break
output.close()
