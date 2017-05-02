import sys, gzip

if len(sys.argv) < 3:
	print 'python input output'
	sys.exit()

output = open(sys.argv[2],"w")

for i in gzip.open(sys.argv[1]):
	if i.startswith("chr"):continue
	x = i.strip("\n").split(",")
	if len(x) < 5:continue ## Comment out, just for testing
	if x[-1] == "NA": continue
	#if float(x[4]) < 0.002:continue
	ratio = float(x[5])

	if ratio >5:
		output.write('chr'+'\t'.join(x)+'\n')
output.close()
