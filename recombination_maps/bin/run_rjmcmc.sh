
## This script runs LDhelmet for whole chromosomes in a single go by processing chunks with a given size. It can be run using the following command:
# parallel -j16 "sh run_whole_chromosome.sh {}" ::: $(ls fasta_files/*)


in_file=$1
prefix=$( echo $in_file | rev | cut -d/ -f1 | rev)
echo $prefix

#exit 0
ldhelmet rjmcmc --num_threads 1 \
		-w 50 \
		-l $2 \
		-p $3 \
		-b 100 \
		-s $in_file \
		--burn_in 100000 \
		-n 1000000 \
		-o processed/$prefix.b100.post \
		-m /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/get_mutation_matrix/non_cpg_sites/rate_matrix.txt \
		-a ancestral/$prefix.anc

ldhelmet post_to_text -m \
		-p 0.025 \
		-p 0.5 \
		-p 0.975 \
		-o processed/$prefix.b100.txt \
		processed/$prefix.b100.post
