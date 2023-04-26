# !/bin/bash
# You need to execute this script in folder with assembly fasta files

mask=$1
curr_dir1="../${mask}/data/${mask}_annotate"
if [[ ! -e ${curr_dir1} ]]
then
        mkdir ${curr_dir1};
fi

curr_dir2="../${mask}/data/orto_rows"
if [[ ! -e ${curr_dir2} ]]
then
        mkdir ${curr_dir2};
fi

curr_dir3="../${mask}/data/for_prokka_fasta"
for f in $(ls ${curr_dir3}/${mask}*.fasta); do
	g=$(echo ${f#${curr_dir3}/});
	#echo ${g%.fasta};
        prokka --outdir ${f%.fasta} cpus 0 --prefix ${g%.fasta} ${f};
done
#
for j in $(ls ${curr_dir3}/${mask}*); do
        if [[ "${j}" == "${j%.*}.gbk" ]]; then
		h=${j%.gbk}
                mv ${curr_dir3}/${h}/${h}.gbk ${curr_dir1}
        fi
done
