### I am not sure whether you need to generate the likelihood lookup table fo reach region individulally, or whether a single one will do. I'll try both and compare...

for i in $1/*
do
python /home/booker/2.Phasing-castaneus/bin/LD_HELMET.py -f $i -p $2 -o $i
done
