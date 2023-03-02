# !/bin/bash

mkdir ../evolution_models
for f in $(ls C_psittaci*.afa); do
	modeltest-ng static -i ${f} -o ../evolution_models/${f%.afa}modeltest;
done


#modeltest-ng static -i C_psittaci_non_uni2_736.afa -o C_psittaci_non_uni2_736.afa_modeltest
