for i in $(ls $1/*.post)
do
prefix=$( echo $i | rev | cut -d/ -f1 | rev)

ldhelmet post_to_text -m \
		-p 0.025 \
		-p 0.5 \
		-p 0.975 \
		-o $prefix.txt \
		$i
done 
