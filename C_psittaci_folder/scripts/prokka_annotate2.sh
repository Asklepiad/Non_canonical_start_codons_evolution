# !/bin/bash
# You need to execute this script in folder with assembly fasta files

mask="C_psittaci"
mkdir ../${mask}_annotate
mkdir ../for_prokka_fasta
mkdir ../orto_rows

for f in $(ls ./${mask}*.fasta); do
        prokka --outdir ./${f%.fasta} --prefix ${f%.fasta} ./${f};
done

for j in $(ls ${mask}*); do
        if [[ "${j}" == "${j%.*}.gbk" ]]; then
                mv ./${j%.*}/${j%.*}.gbk ../${mask}_annotate
        fi
done
