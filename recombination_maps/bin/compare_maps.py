import pandas as pd
import numpy as np
import sys, argparse
import matplotlib.pyplot as plt
import seaborn as sns

def roundto(number, multiple):
	return number+multiple/2 - ((number+multiple/2) % multiple)

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input",
        required = True,
	type = str,
	help = "the processed output of LDhelmet")
parser.add_argument("-b", "--bucket",
        required = False,
	type = int,
	help = "by what number number do you want to group the data? i.e. you will round positions to the nearest multiple of this value ",
	default = 5000)
parser.add_argument("-c", "--chr",
        required = True,
	type = str,
	help = "Give the chromosome that thes segments come from ")
parser.add_argument("-s", "--start",
        required = False,
	type = int,
	help = "Give the starting position of this chunk",
	default = 0)
parser.add_argument("-o", "--output",
        required = True,
	type = str,
	help = "The name of the output CSV")

args = parser.parse_args()


results = pd.read_csv(args.input)
bruns_full = pd.read_csv("/home/booker/mouse_genome/recombination_maps/brunschwig_S1.csv")
bruns_full.columns = ["Chr","Kb","rho_kb"]


bruns =  bruns_full[bruns_full.Chr == args.chr]

bruns["position"] = roundto((bruns.Kb*1000),args.bucket)

results.columns = [ "left_snp","right_snp","mean_rho","p0.025","p0.500","p0.975"]

results["position"] = roundto(((args.start+results.right_snp) + (args.start+results.left_snp))/2,args.bucket)



cast_regions =  {}
for pos, region_data in results.groupby("position"):
	cast_regions[int(pos)] = region_data.mean_rho.median()


bruns_regions =  {}
for pos, bruns_data in bruns.groupby("position"):
	if pos > results.position.max():continue
	bruns_regions[int(pos)] = bruns_data.rho_kb.median()
#pd.DataFrame.from_dict(bruns_regions)


#df = pd.DataFrame([[col1,col2] for col1, d in bruns_regions.items() for col2 in d.items()])
cast =  pd.DataFrame(cast_regions.items(),columns = ["Position (bp)","rho"])
cast["source"] = "cast"
bruns =  pd.DataFrame(bruns_regions.items(),columns = ["Position (bp)","rho"])
bruns["source"] = "bruns"
output = pd.concat([cast,bruns])

for mappy, data in output.groupby("source"):
	print "MEAN"
	print mappy,data.rho.mean()
	print "MEDIAN"
	print mappy,data.rho.median()

output.to_csv(args.output)

