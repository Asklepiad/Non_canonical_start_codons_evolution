# Scripts, their function and options

## Table of content

## Master scripts

## Parts of pipeline

### Parsing_NCBI_1.py
```Parsing_NCBI_1.py``` - creates a list of genbank assemblies with annotations and data about plasmid files (it is important because information about the type of DNA molecules will vanish on the next stage).
### prokka_annotate2.sh
```prokka_annotate2.sh``` - consistently reannotates assemblies. You can read more about the Prokka tool on [its github](https://github.com/tseemann/prokka). We used version 1.14.6.
### First_table_creating3.py
```First_table_creating3.py``` - creates a table about gene features and fasta files with all genes.
### proteinortho_script4.sh
```proteinortho_script4.sh``` - computes ortologous rows and creates a table with information about them. You can read more about the Proteinortho tool on [its gitlub](https://gitlab.com/paulklemm_PHD/proteinortho). We used version 6.1.7.
### Muscle_preparing_5.py
```Muscle_preparing_5.py``` - creates two big summary tables about genes and ortologous rows.
### Statscript.R
```Statscript.R``` - computes different statistics, creates figures and a short report.
### muscle_align6.sh
```muscle_align6.sh``` or ```prank_align6.sh``` aligns sequences. You can read more about the MUSCLE and PRANK aligners respectively on the [MUSCLE v.5 website](https://www.drive5.com/muscle/) and [PRANK website](http://wasabiapp.org/software/prank/). We used version 5.1 and v.170427 of MUSCLE and PRANK respectively.
### ```modeltree_maker7.sh```
```modeltree_maker7.sh``` **(сделать ссылки)** - chooses the best evolution model for further tree creating. You can read more about the modeltest-ng tool on [its github](https://github.com/ddarriba/modeltest). We used version 0.1.7.
### raxml_tree8.sh
```raxml_tree8.sh``` **(сделать ссылки)** - creates phylogenetic trees in Newick format. You can read more about the RAxML-ng tool on [its github](https://github.com/amkozlov/raxml-ng). We used version 1.1.0.
### Ete3_maker_10.py
```Ete3_maker_10.py``` - visualizes phylogenetic trees (*unavailable in server version now*)

## Additional scripts

## Notebooks
