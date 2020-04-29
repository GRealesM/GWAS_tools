#!/bin/bash

shopt -s nullglob
array=(*.gz)

for f in ${array[@]};
do
	FILEBASENAME=$(echo $f | cut -d. -f1); # Take as file base name everything before the first dot
	echo "Working on $FILEBASENAME file"
	zcat $f > temp1
	N_SPACES=$(cat temp1 | head -2 | tail -1 | tr -d -c ' ' | wc -m) # Detect the delimiter by counting the number of spaces on the second line
	if [ $N_SPACES -gt 0 ] 
	   then
#	   echo $N_SPACES
	   sed -e 's/^ *//' -e 's/ *$//' -e 's/ \{1,\}/\t/g' temp1 > "${FILEBASENAME}-fc.tsv" 
	else
#	   echo $N_SPACES
	   sed -e 's/,/\t/g' -e 's/;/\t/g' -e 's/ \{1,\}/_/g' temp1 > "${FILEBASENAME}-fc.tsv"
	fi	
	echo "Field correcting done, now zipping..."
	gzip "${FILEBASENAME}-fc.tsv"
	rm temp1
	mv "${FILEBASENAME}-fc.tsv.gz" ../02-Field_corrected/
done
echo "Done, your files are now at 02-Field_corrected directory"
