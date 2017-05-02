
## This script runs LDhelmet for whole chromosomes in a single go by processing chunks with a given size. It can be run using the following command:
# parallel -j16 "sh run_whole_chromosome.sh {}" ::: $(ls fasta_files/*)


in_file=$1
prefix=$( echo $in_file | rev | cut -d/ -f1 | rev)
echo $prefix
### this doozy takes all the .conf files in the 'processed' directory and tries to run LDhelmet using each one. When it can't get any to work it makes the config files fresh

 
(for i in processed/*conf;do prefix2=$( echo $i | cut -d. -f1 ); sh /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/recombination_maps/bin/run_rjmcmc.sh $in_file $prefix2.lk $prefix2.pade  ; done) || (sh /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/recombination_maps/bin/run_all_ldh_steps.sh $in_file)
