#!/bin/bash

json_path=$1
aligner=$2
po_parameter=$3
mail=$4
threads=$5

if [[ ! -e ../data/jsons ]]
then
        mkdir ../data/jsons;
fi

conda run -n python_start_codons python3 Organisms_parsing02.py "$json_path"

for f in $(ls ../data/jsons); do
	echo "$json_path" "yo-ho-ho!"
        ./folder_creators_server0.sh "$f" "$aligner" "$po_parameter" "$mail" "$threads";
done
