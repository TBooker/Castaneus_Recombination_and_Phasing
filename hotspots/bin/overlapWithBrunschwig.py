
brunschwig ='/home/booker/project/6.Phasing_Castaneus+Recombination_Rates/hotspots/bin/TableS2.csv'
c = 0

for i in open(brunschwig):
	if c==0:
		c=1
		continue
	x= i.strip().split(',')
	print '\t'.join(x)
		
