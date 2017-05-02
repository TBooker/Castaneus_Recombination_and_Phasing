echo "Usage:"
echo "TABIX_RANGES VCf_FILE REFERENCE_GENOME"

# Shell to get the fasta input for LDhelmet...
# Give the name of the tabix file list at the command line as the FIRST input
# Give the name of the VCF file at the command line as the SECOND input
# Give the name of the reference file at the command line as the THIRD input
#mkdir fasta_files
mkdir fasta_files
while read p; do
	echo $p
	SubString=$(echo $p |sed 's/:/_/')
	python /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/ShapeIt2/bin/vcf2fasta.lanner.py -v $2 -s all -i $p -r $3 --sample_ploidy 2 --phase_bool > fasta_files/$SubString
	
done < $1 




