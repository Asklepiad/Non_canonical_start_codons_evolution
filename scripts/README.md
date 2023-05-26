# Scripts, their function and options

## Table of content

[Master scripts](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#master-scripts)

[Parts of pipeline](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#parts-of-pipeline)

[Parsing_NCBI_1.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#parsing_NCBI_1.py)

[prokka_annotate2.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#prokka_annotate2.sh)

[First_table_creating3.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#first_table_creating3.py)

[proteinortho_script4.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#proteinortho_script4.sh)

[Muscle_preparing_5.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#muscle_preparing_5.py)

[Statscript.R](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#statscript.R)

[muscle_align6.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#muscle_align6.sh)

[prank_align6.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#prank_align6.sh)

[modeltree_maker7.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#modeltree_maker7.sh)

[raxml_tree8.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#raxml_tree8.sh)

[Ete3_maker_10.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#ete3_maker_10.py)

[Additional scripts](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#additional-scripts)

[Notebooks](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#notebooks)


https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#master-scripts

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
### modeltree_maker7.sh
```modeltree_maker7.sh``` **(сделать ссылки)** - chooses the best evolution model for further tree creating. You can read more about the modeltest-ng tool on [its github](https://github.com/ddarriba/modeltest). We used version 0.1.7.
### raxml_tree8.sh
```raxml_tree8.sh``` **(сделать ссылки)** - creates phylogenetic trees in Newick format. You can read more about the RAxML-ng tool on [its github](https://github.com/amkozlov/raxml-ng). We used version 1.1.0.
### Ete3_maker_10.py
```Ete3_maker_10.py``` - visualizes phylogenetic trees (*unavailable in server version now*)

## Additional scripts

```check_strains.py```
Script check_strains.py takes a path to alignments as input. First, the function loops through each file with alignments and for each file extracts the most present start codon in the aligned genes. 
Then the script creates a dictionary, where the key is the name of the strain (assembly), and the value is the number of genes whose start codon is minor in relation to the majority in the orthological rows. And next it is created a dictionary, where the key is the name of the strain (assembly), and the value is the list of values: number of genes with gaps at the beginning in the strain, number of all genes in the strain and number of genes whose start codon is minor in relation to the majority in its orthological rows.
After this, the function creates the dataframe, where for each strain, in addition to information about the above three parameters in the dictionary, the percentage of genes with gaps and minor start codons from the total number of genes in the strain is indicated.
We have established a threshold: assemblies whose percentage of genes with gaps at the beginning and genes with minor start codons is higher than the median plus 2 standard deviation are written to a separate list that is printed when the function is executed, and subsequently excluded from the analysis (as they are classified by us as low-quality assemblies). 
Also, when the function is executed, 2 graphs of the distribution of the percentage of genes with gaps at the beginning and the percentage of genes with minor start codons for each strain are drawn.

```download_strains.py```
Script download_strains.py accepts a list with the names of bacteria, downloads genome-wide and scaffold assemblies of bacteria, and creates a dataframe, where for each bacterium its Taxonomy ID in NCBI, the number of genome-wide assemblies, the number of scaffold assemblies and the total number of genome-wide and scaffold assemblies are indicated.
Next, the script divides all bacteria into 3 groups:
1 - bacteria with more than 100 whole genome assemblies in NCBI; 
2 - bacteria with less than 100 whole genome assemblies in NCBI, but more than 100 whole genome plus scaffold assemblies; 
3 - bacteria with less than 100 whole genomes plus scaffold assemblies in NCBI.

Then a dictionary is created, where the key is the name of the bacterium, and the value is the list of assembly IDs for each bacterium from groups 2 (all whole genome assemblies and part of scaffold assemblies up to a total number of assemblies equal to 100) and 3 (all scaffold assemblies in NCBI).
Bacteria from group 1 (with more than 100 whole genome assemblies) are run through the PanACoTA pipeline to filter out related assemblies using the Mash genetic distance. Then, if there are less than 100 strains left after filtering, the script downloads data for all filtered assemblies; if more than 100 strains are left after filtering, the script downloads 100 random filtered assemblies.

As a result the script returns a json file with lists of required amount Taxonomy IDs in NCBI for each bacterium.

## Notebooks
