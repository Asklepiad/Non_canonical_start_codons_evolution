# !/bin/bash


# Conda environments

prokka_check=$(conda env list | grep "prokka_start_codons" | wc -l)
if [[ $prokka_check -eq 0 ]]
then
        conda create -n prokka_start_codons
        conda run -n prokka_start_codons conda install -c bioconda prokka==1.14.6
fi

protortho_check=$(conda env list | grep "proteinortho_start_codons" | wc -l)
if [[ $protortho_check -eq 0 ]]
then
        conda create -n proteinortho_start_codons
        conda run -n proteinortho_start_codons conda install -c bioconda proteinortho==6.1.7
fi

modeltest_check=$(conda env list | grep "modeltest_start_codons" | wc -l)
if [[ $modeltest_check -eq 0 ]]
then
        conda create -n modeltest_start_codons
        conda run -n modeltest_start_codons conda install -c bioconda modeltest-ng==0.1.7
fi

raxml_check=$(conda env list | grep "raxml_start_codons" | wc -l)
if [[ $raxml_check -eq 0 ]]
then
        conda create -n raxml_start_codons
        conda run -n raxml_start_codons conda install -c bioconda raxml-ng==1.1.0
fi

python_check=$(conda env list | grep "python_start_codons" | wc -l)
if [[ $python_check -eq 0 ]]
then
        conda create -n python_start_codons python==3.9.13
        conda run -n python_start_codons pip3 install -r requirements.txt
fi

prank_check=$(conda env list | grep "prank_start_codons" | wc -l)
if [[ $prank_check -eq 0 ]]
then
        conda create -n prank_start_codons
        conda run -n prank_start_codons conda install -c bioconda prank==170427
fi

R_check=$(conda env list | grep "R_start_codons" | wc -l)
if [[ $R_check -eq 0 ]]
then
        conda create -n R_start_codons
        conda run -n R_start_codons conda install -c conda-forge r-base==4.2.2
fi

muscle_check=$(conda env list | grep "muscle_start_codons" | wc -l)
if [[ $muscle_check -eq 0 ]]
then
        conda create -n muscle_start_codons
        conda run -n muscle_start_codons conda install -c bioconda muscle==5.1
fi

