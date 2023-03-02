# !/bin/bash

directory="./multialignments/"
for f in $(ls ${directory}*.fasta); do
	length=$(cat $f | grep -A 1 ">" - | head -2 | wc -c);
	if [[ $length -lt 1000 ]]
	then
		muscle -align $f -output ${f%.fasta}.afa;
	else
		muscle -super5 $f -output ${f%.fasta}.afa;
	fi
done
