import glob, sys
if len(sys.argv) < 3:
	print "You give this script the path to a directory that contains a bunch of post_to_text output from LDhelmet and makes a concatenated recombination map\n\tpython\tmake_chr_map.py\tpath_to_dir\toutput_name"
output = open(sys.argv[2],"w")

for i in glob.glob(sys.argv[1]+"/*b10.txt"):
#	print i
	start = i.split("/")[-1].split(".")[0].split("_")[1].split("-")[0]
	count = 0

#	print len(open(i).readlines())
	for l in open(i).readlines()[203:-200]:
		
#		break
		count +=1
		if l.startswith("#") or l.startswith("version"):continue
		x = l.strip("\n").split(" ")
		
		l_snp = int(x[0]) + int(start)
		r_snp = int(x[1]) + int(start)
		outline = ",".join(map(str,[l_snp,r_snp]+x[2:]))
		output.write(outline+"\n")
	print i
	print count
	#break

output.close()
