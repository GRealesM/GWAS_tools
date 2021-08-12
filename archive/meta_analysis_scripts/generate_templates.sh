#!/bin/bash
# Generating input files for meta-analysis

traits=($(ls *.gz | cut -d_ -f1 | sort | uniq)) # I know that piping ls output is BAD, but I did for for practicity

for index in ${!traits[@]};
do
	filenames=($(ls *.tsv.gz | grep "${traits[index]}" | tr "\n" " "))
	inputfile1=${filenames[0]}
	inputfile2=${filenames[1]}
	sed METAL_template.txt -e 's/inputfile1.txt/'"$inputfile1"'/' -e 's/inputfile2.txt/'"$inputfile2"'/' -e 's/_inputtrait/_'"${traits[index]}"'/' > METAL_input_${traits[index]}
done
