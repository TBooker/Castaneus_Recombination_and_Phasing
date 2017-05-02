### Quick python script to convert LDhelmet output to a csv file

import sys

output = open(sys.argv[1].replace(".txt",".csv"),"w")

for i in open(sys.argv[1]):
	if i.startswith("# parameters"):
		continue
	elif i.startswith("version"):
		continue
	elif i.startswith("# left"):
		output.write(",".join(i.split(" ")[1:]))
	else:
		output.write(",".join(i.split(" ")))
output.close()