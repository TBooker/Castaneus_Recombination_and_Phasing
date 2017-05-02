import subprocess, argparse
import pandas as pd

def main():
	parser = argparse.ArgumentParser(description="Takes a bed file, shuffles it and sees how much overlap there is between the shuffled .bed and a reference .bed . Uses bedtools to acheive this")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input", 
		type =str, 
		help = "The input .bed file that you want to shuffle and compare")

	parser.add_argument("-c","--comparison", 
		required = True,
		dest = "comparison", 
		type =str, 
		help = "The .bed file that you want to compare to")

	parser.add_argument("-g","--genome_file", 
		required = True,
		dest = "genome_file", 
		type =str, 
		help = "The genome file used by bedtools")

	parser.add_argument("-e","--exclude", 
		required = False,
		dest = "exclude", 
		type =str, 
		help = "Regions of the genome you want to exclude from being shuffled (e.g. the first 3e6 bases on most of mm9 chromosomes",
		default = '')

	parser.add_argument("-r","--reps", 
		required = True,
		dest = "reps", 
		type =int, 
		help = "reps")

	parser.add_argument("-o","--output", 
		required = True,
		dest = "output", 
		type =str, 
		help = "Name the output file (will be a csv file)")

	args = parser.parse_args()
	
	results = []	

	
	exclude_line = ''
	for rep in range(args.reps):

		shuffle_command = 'bedtools shuffle -chrom -i ' + args.input + ' -g ' + args.genome_file + exclude_line + '> shuffled.'+str(rep)+'.temp.bed'

		shuffle = subprocess.Popen(shuffle_command, stdout=subprocess.PIPE, stderr=None, shell=True)

		shuffle.communicate()

		comparison_command = 'bedtools intersect -wa -a ' + args.comparison + ' -b shuffled.'+str(rep)+'.temp.bed | wc -l'

		comparison = subprocess.Popen(comparison_command, stdout=subprocess.PIPE, stderr=None, shell=True)

		overlap = comparison.communicate()[0].strip()

		subprocess.Popen('rm shuffled.'+str(rep)+'.temp.bed', stdout=subprocess.PIPE, stderr=None, shell=True).communicate()

		results.append( [rep, overlap] )

	pd.DataFrame(results, columns = ['rep', 'overlap']).to_csv(args.output, index =False)

if '__name__':
	main()
