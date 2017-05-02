import sys
#chromosome_lengths = {1:197195432,10:130694993,11:122082543,12:120129022,13:120421639,14:124902244,15:104043685,16:98207768,17:94987271,18:90702639,19:61431566,2:182113224,3:160039680,4:156508116,5:151834684,6:149736546,7:145441459,8:129401213,9:124595110}
#chromo_rates = {1:0.61e-8, 2:0.54e-8, 3:0.44e-8, 4:0.54e-8, 5:0.56e-8, 6:0.45e-8, 7:0.52e-8, 8:0.61e-8, 9:0.58e-8, 10:0.61e-8, 11:0.68e-8 ,12:0.51e-8, 13:0.52e-8, 14:0.64e-8, 15:0.63e-8, 16:0.54e-8, 17:0.57e-8, 18:0.46e-8, 19:1.04e-8}

rho_sum = 0 
bases = 0
previous_pos = 0
for i in open(sys.argv[1]):
	x = i.split("\t")
#	print x
	rho = float(x[6])
	pos = int(x[3])

	bases+= (pos-previous_pos)
	rho_sum +=(pos-previous_pos)*rho		
	previous_pos = pos	
bases+=(1000000-previous_pos)
rho_sum +=(1000000-previous_pos)*rho

print bases
print rho_sum
print (rho_sum/bases)/4

#outlist_c = map(str,["castaneus",c,bases,rho_sum,((rho_sum/bases)/chromo_rates[c])/4])
#output.write(",".join(outlist_c)+"\n")



