# Scripts, their function and options

## Table of content

[Master scripts](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#master-scripts)

[Parts of pipeline](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#parts-of-pipeline)

[Parsing_NCBI_1.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#parsing_NCBI_1py)

[prokka_annotate2.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#prokka_annotate2sh)

[First_table_creating3.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#first_table_creating3py)

[proteinortho_script4.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#proteinortho_script4sh)

[Muscle_preparing_5.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#muscle_preparing_5py)

[Statscript.R](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#statscriptR)

[muscle_align6.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#muscle_align6sh)

[prank_align6.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#prank_align6sh)

[modeltree_maker7.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#modeltree_maker7sh)

[raxml_tree8.sh](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#raxml_tree8sh)

[Ete3_maker_10.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#ete3_maker_10py)

[Additional scripts](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#additional-scripts)

[check_strains.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#check_strainspy)

[download_strains.py](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#download_strainspy)

[Notebooks](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#notebooks)

## Master scripts

## Parts of pipeline

#### ```Parsing_NCBI_1.py```
> Short description. 

Сreates a list of genbank assemblies with annotations and data about plasmid files (it is important because information about the type of DNA molecules will vanish on the next stage).

> Input. 

1. Path to json file, which named equally to organism of interest name's long form.
2. Your e-mail address.

> Output. 

1. Content of "for_prokka_fasta" folder (fasta files for given assemblies).
2. Two json files with the data about which file contents plasmid fasta and which has chromosomal one.

> Detailed description. 

Script finds all links to scaffold and complete genome assemblies, checks them for unicity, collected all annotations for assemblies (gbk files). On the next stage it saved sequences into fasta files for reannotating on the further stage. In addition script creates two json files, which content data about DNA source of each fasta file: if it is a chromosome or plasmid.

#### ```prokka_annotate2.sh```
> Short description. 

Consistently reannotates assemblies. 

> Input. 

1. Short name of the bacteria (first letter of genus name, underscore, specie's name or first nine letters of specie's name, if it longer). The shortening of specie's name implemented because long second name lead to the errors.

> Output. 

Annotations (with gbk extensions) in `<organism_annotate>` directory

> Version and options.

We used version 1.14.6. 

Options:

1. ```cpus 0``` means prokka will use all available threads. If you want to change behavior of prokka, you may create new bash variable in prokka script (and mention it in  `folder_creators_server0.sh`) or explicitly write number of cpus in prokka script.

> Detailed description.

You can read more about the Prokka tool on [its github](https://github.com/tseemann/prokka).

#### ```First_table_creating3.py```
> Short description. 

Creates a table about gene features and fasta files with all genes.

> Input. 

1. Short name of the bacteria

> Output.

1. Table with data about each gene: id, assembly, nucleotide sequence, aminoacid sequence, type od DNA molecule (plasmid or chromosome), start-codone, protein product of gene (if known), COG (functionality category of gene), coordinates of start and end of the gene on chromosome, strand, length of gene, special variable for filtering uncomplete sequences, accssory status to each of 25 COG categories.
2. Fasta files of each selected assembly for further analysis.

> Detailed description.

Script pasrses assemblies annotations (gbk files from prokka) and creates big table, where samples are genes and features are some characteristics (described in point one of this scripts's `output`). Table is in csv format. Table's detailed desxription you can find in `../data` folder's readme file. For assigning functional category(ies) for each gene script uses COG tables from `../data/bigcog/` folder.
In addition script wtites fasta files for each assembly. It needs for further ortologous groups construction.


#### ```proteinortho_script4.sh```
> Short description. 

Computes ortologous rows and creates a table with information about them. You can read more about the Proteinortho tool on [its gitlub](https://gitlab.com/paulklemm_PHD/proteinortho). We used version 6.1.7.

> Input. 

> Output. 

> Detailed description.


#### ```Muscle_preparing_5.py```
> Short description. 

Creates two big summary tables about genes and ortologous rows.

> Input. 

> Output. 

> Detailed description.


#### ```Statscript.R```
> Short description. 

Computes different statistics, creates figures and a short report.

> Input. 

> Output. 

> Detailed description.


#### ```muscle_align6.sh``` or ```prank_align6.sh```
> Short description. 

Aligns sequences. You can read more about the MUSCLE and PRANK aligners respectively on the [MUSCLE v.5 website](https://www.drive5.com/muscle/) and [PRANK website](http://wasabiapp.org/software/prank/). We used version 5.1 and v.170427 of MUSCLE and PRANK respectively.

> Input. 

> Output. 

> Detailed description.


#### ```modeltree_maker7.sh```
> Short description. 

Chooses the best evolution model for further tree creating. You can read more about the modeltest-ng tool on [its github](https://github.com/ddarriba/modeltest). We used version 0.1.7.

> Input. 

> Output. 

> Detailed description.


#### ```raxml_tree8.sh```
> Short description. 

Creates phylogenetic trees in Newick format. You can read more about the RAxML-ng tool on [its github](https://github.com/amkozlov/raxml-ng). We used version 1.1.0.

> Input. 

> Output. 

> Detailed description.


#### ```Ete3_maker_10.py```
> Short description. 

Visualizes phylogenetic trees (*unavailable in server version now*)

> Input. 

> Output. 

> Detailed description.


## Additional scripts

#### ```check_strains.py```
> Short description. 

> Input. 

> Output. 

> Detailed description.

Script check_strains.py takes a path to alignments as input. First, the function loops through each file with alignments and for each file extracts the most present start codon in the aligned genes. 
Then the script creates a dictionary, where the key is the name of the strain (assembly), and the value is the number of genes whose start codon is minor in relation to the majority in the orthological rows. And next it is created a dictionary, where the key is the name of the strain (assembly), and the value is the list of values: number of genes with gaps at the beginning in the strain, number of all genes in the strain and number of genes whose start codon is minor in relation to the majority in its orthological rows.
After this, the function creates the dataframe, where for each strain, in addition to information about the above three parameters in the dictionary, the percentage of genes with gaps and minor start codons from the total number of genes in the strain is indicated.
We have established a threshold: assemblies whose percentage of genes with gaps at the beginning and genes with minor start codons is higher than the median plus 2 standard deviation are written to a separate list that is printed when the function is executed, and subsequently excluded from the analysis (as they are classified by us as low-quality assemblies). 
Also, when the function is executed, 2 graphs of the distribution of the percentage of genes with gaps at the beginning and the percentage of genes with minor start codons for each strain are drawn.

#### ```download_strains.py```
> Short description. 

> Input. 

> Output. 

> Detailed description.

Script download_strains.py accepts a list with the names of bacteria, downloads genome-wide and scaffold assemblies of bacteria, and creates a dataframe, where for each bacterium its Taxonomy ID in NCBI, the number of genome-wide assemblies, the number of scaffold assemblies and the total number of genome-wide and scaffold assemblies are indicated.
Next, the script divides all bacteria into 3 groups:
1 - bacteria with more than 100 whole genome assemblies in NCBI; 
2 - bacteria with less than 100 whole genome assemblies in NCBI, but more than 100 whole genome plus scaffold assemblies; 
3 - bacteria with less than 100 whole genomes plus scaffold assemblies in NCBI.

Then a dictionary is created, where the key is the name of the bacterium, and the value is the list of assembly IDs for each bacterium from groups 2 (all whole genome assemblies and part of scaffold assemblies up to a total number of assemblies equal to 100) and 3 (all scaffold assemblies in NCBI).
Bacteria from group 1 (with more than 100 whole genome assemblies) are run through the PanACoTA pipeline to filter out related assemblies using the Mash genetic distance. Then, if there are less than 100 strains left after filtering, the script downloads data for all filtered assemblies; if more than 100 strains are left after filtering, the script downloads 100 random filtered assemblies.

As a result the script returns a json file with lists of required amount Taxonomy IDs in NCBI for each bacterium.

## Notebooks
