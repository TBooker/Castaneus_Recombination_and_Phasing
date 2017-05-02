cd chr$1/

#mkdir fasta_files
#mkdir ancestral
#mkdir processed
#mkdir map

#python /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/ShapeIt2/bin/anc_allele.py -v /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/ShapeIt2/per_chromosome/chr$1/chr$1.phased.vcf.gz -m  /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/get_mutation_matrix/non_cpg_sites/rate_matrix.txt -f chr$1.tabixRanges.txt -c $1 -o ancestral/

#sh ~/project/6.Phasing_Castaneus+Recombination_Rates/ShapeIt2/bin/get_fastas.sh chr$1.tabixRanges.txt ~/project/6.Phasing_Castaneus+Recombination_Rates/ShapeIt2/per_chromosome/chr$1/chr$1.phased.vcf.gz /home/booker/mouse_genome/mouse_genome_build_mm9/MAF_castaneus_swapped/new.chr$1.fa

parallel -j15 "sh /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/recombination_maps/bin/run_whole_chromosome_old.sh {}" ::: $(ls fasta_files/*)

python /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/recombination_maps/bin/make_chr_map.py processed map/chr$1.b10.map

sort -nk1 map/chr$1.b10.map > map/chr$1.b10.sorted.map

python /home/booker/mouse_genome/bin/condense_ldhelmet_output.py map/chr$1.b10.sorted.map chr$i map/chr$1.b10.condensed.map

rm map/chr$1.map

cd ../
