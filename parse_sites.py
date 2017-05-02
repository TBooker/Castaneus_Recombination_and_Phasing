import sys
## quick script to parse te  sites file output of HapSeq into a format accepted by LDhelmet
if len(sys.argv) < 2:
	print "Need to give the sites file name at the command line"
	sys.exit()

new_file = open(sys.argv[1][:-4]+"parsed.sites","w")

with open(sys.argv[1]) as FILE:
	for line in FILE:
		x = line.split("\t")
		new_file.write(x[0]+"\n")

