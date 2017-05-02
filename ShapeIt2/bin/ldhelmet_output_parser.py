import glob, sys, os, argparse
from multiprocessing import Pool
from tom import brace
parser = argparse.ArgumentParser(description="This script will combine and simplify LDhelmet output for a single chromosome. For example, severl places may have large regions with constant recombination rate, these do not all need a separate line so merge these.")

parser.add_argument("-d","--directory", 
		required =True,
		dest = "directory", 
		type =str, 
		help = "The directory that contains all the post files from LDhelmet")
parser.add_argument("-o","--output", 
		required =True,
		dest = "output", 
		type =str, 
		help = "The name you want to give to the resulting file")
parser.add_argument("--new_dir", 
		required =False,
		dest = "new_dir", 
		type =str, 
		help = "The name you want to give to the folder conataining the LDhelmet post-to-text output",
		default="ldhelmet_processed")
parser.add_argument("-p","--procs", 
		required =False,
		dest = "procs", 
		type =int, 
		help = "How many threads are you happy to split th process onto?",
		default=1)
parser.add_argument("--skip", 
		required =False,
		dest = "skip", 
		action = "store_true", 
		help = "Have you already got the files you need? For testing purposes")
####  FUNCTIONS...

def process_ldhelmet(file_name):
	name = file_name.split("/")[-1].split(".")[0]+".rho.txt"
	os.system("ldhelmet post_to_text -m -o "+name+" " + file_name)
## Could add different things to this file

def process_output(file_name):
	start = int(file_name.split("/")[-1].strip(".rho.txt").split("_")[1].split("-")[0])
	regions =[]
	for i in open(file_name,"r"):
		if i.startswith("#"):continue
		elif i.startswith("version"):continue
		x = i.strip("\n").split(" ")
		pos = int(x[0])+start
		rho = float(x[2])
		if len(regions) == 0:
			regions.append([pos,rho])
		elif regions[-1][1] == rho: continue
		else:
			regions.append([pos,rho])
	return regions			
		

####

args = parser.parse_args()

if not args.skip:
	os.system("mkdir " + args.new_dir)
	file_list = glob.glob(args.directory + "/*.post")
	if len(file_list)==0:
		print "The directory you specified does not contain any .post files"
		sys.exit()
	
	print "there are",len(file_list),"files to process"
	## Is there the need to multiprocess this??
	if args.procs == 1:
		[process_ldhelmet(i) for i in file_list]
	elif args.procs > 1:
		p = Pool(args.procs)
		p.map(process_ldhelmet,file_list)
	
	os.system("mv *.rho.txt " + args.new_dir +"/")
#########################################
########	Act II		#########
#########################################

file_list = glob.glob(args.new_dir+"/*.rho.txt")
output = open(args.output,"w")
for i in file_list:
	vals = process_output(i)
	for j in vals:
		output.write(",".join(map(str,j))+"\n")













