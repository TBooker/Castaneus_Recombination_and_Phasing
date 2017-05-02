#
### Before running this script make sure that there is a directory called 'processed' in the directory
#
#
## This script runs LDhelmet for whole chromosomes in a single go by processing chunks with a given size. It can be run using the following command:
# parallel -j16 "sh run_whole_chromosome.sh {}" ::: $(ls fasta_files/*)


in_file=$1
prefix=$( echo $in_file | rev | cut -d/ -f1 | rev)
echo $prefix

#exit 0
#ldhelmet find_confs --num_threads 1 \
#		-w 50 \
#		-o processed/$prefix.conf \
#		$in_file
#ldhelmet pade --num_threads 1 \
#		-c processed/$prefix.conf \
#		-t 0.01 \
#		-x 11 \
#		-o processed/$prefix.pade
#ldhelmet table_gen --num_threads 1 \
#		-c processed/$prefix.conf \
#		-t 0.01 \
#		-r 0.0 0.1 10.0 1.0 100.0 \
#		-o processed/$prefix.lk
ldhelmet rjmcmc --num_threads 1 \
		-w 50 \
		-l processed/$prefix.lk \
		-p processed/$prefix.pade \
		-b 100 \
		-s $in_file \
		--burn_in 100000 \
		-n 1000000 \
		-o processedNoAncAllele/$prefix.b100_NoAncAllele.post \
		-m /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/get_mutation_matrix/non_cpg_sites/rate_matrix.txt \
#		-a ancestral/$prefix.anc

ldhelmet post_to_text -m \
		-p 0.025 \
		-p 0.5 \
		-p 0.975 \
		-o processedNoAncAllele/$prefix.b100_NoAncAllele.txt \
		processedNoAncAllele/$prefix.b100_NoAncAllele.post
