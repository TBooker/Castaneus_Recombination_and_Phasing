for i in $(seq 6 19)
do
echo $i
python /home/booker/project/6.Phasing_Castaneus+Recombination_Rates/ShapeIt2/bin/danVCF_2_oxfordVCF.py -v /home/booker/mouse_genome/vcf_files/regular/chr$i.calls.vcf.gz -c chr$i --HWE --minDP 10 --minGQ 15 --minQUAL 30 -o chr$i
bgzip chr$i.ox.vcf
done
