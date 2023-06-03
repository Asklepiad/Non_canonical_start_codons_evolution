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

There are two modes of pipeline working:

- You can start pipeline from parsing "big" json file, which is output of `Download_strains_03.py`. Then you will serially run the pipeline on each bacteria from json. You need to use `json.sh` script for it.
**пример команды**
- On another way, you may run pipeline on one given bacteria. It may be realized by running `folder_creators_server0.sh`.
**пример команды**

![Big and small pipelines](./illustrations/big_small_pipeline.png)

The more detailed description will be below.

#### `json.sh`

> Input. 

1. Path to *many bacterias* json file.
2. Aligner (muscle or prank in current realisation)
3. Identity percent of best blast hits in proteinortho.
4. Your e-mail
5. Number of threads for parallelising some of tools.

> Output.

No explicit output. Runs `Organisms_parsing02.py` and `folder_creators_server0.sh` (the last one described in the section below)
In fact results of running are equivalent to `folder_creators_server0.sh`, but for more than one specie.

#### `folder_creators_server0.sh`

> Input. 

1. Path to *one bacteria* json file.
2. Aligner (muscle or prank in current realisation)
3. Identity percent of best blast hits in proteinortho.
4. Your e-mail
5. Number of threads for parallelising some of tools.

> Output.

Directory, which named equivalently to shortened bacteria name. Contents two subdirectories:
1. THe data folder contents csvs, alignments, different fasta files, temporary jsons, results of tools' computings etc
2. The figures folder contents plots and phylogenetic trees (in separate directory).


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

Annotations (with gbk extensions) in `../<organism>/data/<organism>_annotate` directory.

> Version and options.

We used version 1.14.6. 

Options:

1. ```cpus 0``` means prokka will use all available threads. If you want to change behavior of prokka, you may create new bash variable in prokka script (and mention it in  `folder_creators_server0.sh`) or explicitly write number of cpus in prokka script.

> Detailed description.

Creates the folders `../<organism>/data/<organism>_annotate/` and `../<organism>/data/for_prokka_fasta/`
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

Script pasrses assemblies annotations (gbk files from prokka) and creates big table, where samples are genes and features are some characteristics (described in point one of this scripts's `output`). Table is in csv format. Table's detailed desxription you can find in `../output_example/data` folder's readme file. For assigning functional category(ies) for each gene script uses COG tables from `../data/bigcog/` folder.
In addition script wtites fasta files for each assembly in the `../<organism>/data/orto_rows/` directory. It needs for further ortologous groups construction.


#### ```proteinortho_script4.sh```
> Short description. 

Computes ortologous rows (OR) and creates a table with information about them. 

> Input. 

1. The short name of the bacteria.
2. Identity percent of best blast hits (our default value is 75, which is a standart, when you work with strains from one specie).

> Output. 

1. Files `S_ruber.proteinortho-graph`, `S_ruber.proteinortho-graph.summary`, `S_ruber.proteinortho.html`, `S_ruber.blast-graph` and `S_ruber.proteinortho.tsv`. We are interested in the last one. More detailed description about each of them you can find in readme of the `../output_example/data` folder.

> Options.

1. ```--debug``` gives detailed information for bug tracking
2. ```--selfblast``` detects paralogs without orthologs
3. ```--singles``` prints singletones into the output file
4. ```--identity``` percent of best blast hits

> Detailed description.

Creates the folder `../<organism>/data/multialignments/`
You can read more about the Proteinortho tool on [its gitlub](https://gitlab.com/paulklemm_PHD/proteinortho). We used version 6.1.7.

#### ```Muscle_preparing_5.py```
> Short description. 

Creates two big summary tables about genes and ORs. Writes fasta files from ORs with different start-codons for further analysis.

> Input. 

1. Path to json file, which named equally to organism of interest name's long form.

> Output. 

1. Fasta files for alignment in `../<organism>/data/multialignments`
2. Table with data about each gene (upgraded version of `First_table_creating3` output).
3. Table with data about each OR.
4. Shortened version of the previous one for one kind of analysis.

> Detailed description.
Script filters ORs with paralogs, assigns OR number for each gene, collects data about start codons proportion in each OR. On the next stage script evaluates uniformity of the OR. It considers uniform, only if all start codons in a row are same. Then script writes All non-uniform rows sequences into the `../<organism>/data/multialignments` folder. 

At the end script combined three tables with data:
- Upgraded version of `First.table.csv` with extended data about each non-paralogous gene.
- Table with data about each OR and its shortened version, from which was excluded information about source of gene in it (number of assembly).

Script has some primitive time-logging, which was explored, when we imporoved the effectiveness of code (first version works *378 minutes and 14 seconds* on one of the bacteria before upgrading and *58 seconds* after upgrading).

#### ```Statscript.R```
> Short description. 

Computes different statistics, creates figures and a short report.

> Input. 

1. Short name of the bacteria

> Output. 

1. `<organism>_statistical_short_report.txt` Short report with statistical data about bacteria: number of genes and ORs, distributions of SCs, proportion of each SCs, SC-cog relationships, data about pangenome fractions distribution for each type of SCs, uniformity of the ORs.
2. `<organism>_for_common_boxplot.csv` Csv table with data of gene number with each SC in each assembly. On the next stage of analysis these table from all species of interest will be used in common boxplot creation.
3. `<organism>_noncanonic_products.csv` Non-canonical genes' products' list.
4. `<organism>_C_psittaci_COG_precomputed.csv` Csv table with data about genes' COGs distribution.
5. `<organism>_C_psittaci_COG_start.csv` Csv table with proporions of each COG in all SCs with confident intervals (CI).
6. `<organism>_cogs.csv` Csv table with number of each COG-SC pair (except genes with unknown function).
7. `<organism>_cog_stat_per_sc.csv` Csv table with data about number and percent of SCs.
8. `<organism>_cog_sc_mwu.csv` Csv table with information about statistical difference between SCs realtive proportions for given COG.
9. Figures and plots.

> Detailed description.

Script used information from tables about genes (`summary_rows_prokka.csv`) and ORs (`start_codons2_prokka.csv`) in the folder of the appropriate organism.

On the next stage script computes pangenome statistics (proportion of core, shell and cloud), connections between SCs, uniformity of the ORs and pangenome fractions.

Then script creates gene frequency spectrum, draws boxplots and violinplots. In addition script creates errorbars for relationship between SCs and pangenome fractions, and between SCs and COGs. Script also draws barplots, boxplots, spectrs and other visualisations of the data.

After this computations script evaluates proportions of COGs with different grouping variables, frequencies of ortologous rows with given number of genes with given SC. Then script computes if SCs relative proportion in one COG is different in different types of SCs.

And as a final accord, script creates short statistical report about an organism.

#### ```muscle_align6.sh``` or ```prank_align6.sh```
> Short description. 

Aligns sequences. We used version 5.1 and v.170427 of MUSCLE and PRANK respectively.

> Input. 

1. The short name of the bacteria. (for both MUSCLE and PRANK)
2. Number of threads (for MUSCLE only). Default value from sbatch script is 24.

> Output. 

1. Alignments in fasta format in the `../<organism>/data/multialignments/` folder.

> Detailed description. 

Creates the folder `../<organism>/data/evolution_models/`
You can read more about the MUSCLE and PRANK aligners respectively on the [MUSCLE v.5 website](https://www.drive5.com/muscle/) and [PRANK website](http://wasabiapp.org/software/prank/).

#### ```modeltree_maker7.sh```
> Short description. 

Chooses the best evolution model for further tree creating.  We used version 0.1.7.

> Input. 

1. The short name of the bacteria.
2. Aligner used om the previous step.

> Output. 

1. Files in the ```../<organism>/data/evolution_models/``` folder: tree file with phylogenetic tree in Newick format, out data about model, log file and ckp (checkpoint) file.

> Detailed description.

You can read more about the modeltest-ng tool on [its github](https://github.com/ddarriba/modeltest).

#### ```raxml_tree8.sh```
> Short description. 

Creates phylogenetic trees in Newick format. We used version 1.1.0.

> Input. 

1. The short name of the bacteria.
2. Aligner used om the previous step.
3. Number of threads. Default value from sbatch script is 24.

> Output. 

1. Files: bestmodel (optimized model parameters for the best-scoring ML tree), startree (starting tree(s) used for maximal likelihood (ML) inference), mltrees(ML trees for each starting tree), besttree (best-scoring ML tree), log file and rba (compressed alignment).

> Detailed description.

Creates two new directories: ```../<organism>/data/raxmlng_trees``` and ```../<organism>/data/configs```.
You can read more about the RAxML-ng tool on [its github](https://github.com/amkozlov/raxml-ng).

#### ```Ete3_maker_10.py```
> Short description.

Visualizes phylogenetic trees (*unavailable in server version now*)

> Input. 

1. The short name of the bacteria.

> Output. 

Phylogenetic tree's image is in the ```../<organism>/figures/trees/``` folder

> Detailed description.

Draws phylogenetic trees. Each ATG-related genes has red color, GTG-related - has green color and TTG-realted - has blue color.

## Additional scripts

#### ```Check_strains_04.py```
> Short description. 

The script analyzes all alignments and displays lists of assemblies in which a large number of genes, when aligned within orthologous rows, had gaps at the beginning and had a minor start codon in relation to the majority.

> Input. 

1. Path to alignments.

> Output. 

2 lists: 
1. A list with the names of assemblies in which a large number of genes, when aligned within orthologous rows, had gaps at the beginning;
2. A list with the names of assemblies in which a large number of genes, when aligned within orthologous rows, had a minor start codon in relation to the majority.

> Detailed description.

Script check_strains.py takes a path to alignments as input. First, the function loops through each file with alignments and for each file extracts the most present start codon in the aligned genes. 
Then the script creates a dictionary, where the key is the name of the strain (assembly), and the value is the number of genes whose start codon is minor in relation to the majority in the orthological rows. And next it is created a dictionary, where the key is the name of the strain (assembly), and the value is the list of values: number of genes with gaps at the beginning in the strain, number of all genes in the strain and number of genes whose start codon is minor in relation to the majority in its orthological rows.
After this, the function creates the dataframe, where for each strain, in addition to information about the above three parameters in the dictionary, the percentage of genes with gaps and minor start codons from the total number of genes in the strain is indicated.
We have established a threshold: assemblies whose percentage of genes with gaps at the beginning and genes with minor start codons is higher than the median plus 2 standard deviation are written to a separate list that is printed when the function is executed, and subsequently excluded from the analysis (as they are classified by us as low-quality assemblies). 
Also, when the function is executed, 2 graphs of the distribution of the percentage of genes with gaps at the beginning and the percentage of genes with minor start codons for each strain are drawn.

#### ```Download_strains_03.py```
> Short description. 

The script checks the number of available complete genome and scaffold assemblies for the input organisms, filters the assemblies if necessary using the PanACoTA pipeline, and saves the json file, where each bacterium corresponds to a list of assembly IDs.

> Input. 

1. List of bacteria.
2. Your mail in NCBI.

> Output. 

1. Json file ("id_lists.json") with lists of assemblies ID's for each bacterium.

> Detailed description.

Script download_strains.py accepts a list with the names of bacteria, downloads genome-wide and scaffold assemblies of bacteria, and creates a dataframe, where for each bacterium its Taxonomy ID in NCBI, the number of genome-wide assemblies, the number of scaffold assemblies and the total number of genome-wide and scaffold assemblies are indicated.
Next, the script divides all bacteria into 3 groups:
1 - bacteria with more than 100 whole genome assemblies in NCBI; 
2 - bacteria with less than 100 whole genome assemblies in NCBI, but more than 100 whole genome plus scaffold assemblies; 
3 - bacteria with less than 100 whole genomes plus scaffold assemblies in NCBI.

Then a dictionary is created, where the key is the name of the bacterium, and the value is the list of assembly IDs for each bacterium from groups 2 (all whole genome assemblies and part of scaffold assemblies up to a total number of assemblies equal to 100) and 3 (all scaffold assemblies in NCBI).
Bacteria from group 1 (with more than 100 whole genome assemblies) are run through the PanACoTA pipeline to filter out related assemblies using the Mash genetic distance. Then, if there are less than 100 strains left after filtering, the script downloads data for all filtered assemblies; if more than 100 strains are left after filtering, the script downloads 100 random filtered assemblies.

As a result the script returns a json file with lists of required amount Taxonomy IDs in NCBI for each bacterium.


#### ```Organisms_parsing02.py```

> Description. 

Script parsed json file with data about each organism of interest and its links to assemblies in NCBI. The script is half-additional, because it is a necessary part of `json.sh` pipeline, but unnecessary of `folder_creators`.

> Input. 

1. Path to json file with assembly links for each bacteria.

> Output. 

1. Creates `jsons` folder in the `../data/` directory and creates single json files for each bacteria from input json in it.


#### ```environments_script_01.sh```

> Description. 

Script creates environments on server or local machine with conda. It is helpful for reproducibility, because script downloads exactly those tools version, which were used in original investigation.


#### ```local_statscript.R```

Script is basically equal to server version [`Statscript.R`](https://github.com/Asklepiad/BI_project_2022/tree/main/scripts#statscriptr). The difference is script created for running in RStudio and has another way for installing dependecies. In addition, it doesn't takes variables a-la argparse mode - the name of bacteria needs to be explicitely passed into the text of script.

#### ```common_statistics.R```

> Description. 

Computes statistics by all dataset of bacterias. Draws barplots and boxplots which rebuts SCs distribution between bacterial species.

> Input. 

1. Path to folder with data about all organisms bolxplots.
2. Path to folder with shortened data about ortologous rows of all organisms.
3. List of bacterias' name.

> Output. 

1. `barplot_common_2scs.png` - barplot with proportions of GTG and TTG SCs. Bars grouped by phylogeny.
2. `barplot_common_2scs_size.png"` - barplot with proportions of GTG and TTG SCs. Bars grouped by genome size.
3. `boxplot_common_2scs.png` - boxplot with data about number of different SCs in assemblies of different organisms.

## Notebooks

Some parts of pipeline may be runned from jupyter-notebooks. Functionality is equal to corresponding scripts, except there is no argparse functionality in notebooks (user need to print variables manually) and there are slightly another ways used. Detailed information about functionality, input and output may be found in the appropriate sections for py scripts [```Organisms_parsing02.ipynb```](), [```Parsing_NCBI_1.ipynb```](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#parsing_NCBI_1py), [```First_table_creating3.ipynb```](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#first_table_creating3py), [```Muscle_preparing_5.ipynb```](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#muscle_preparing_5py), [```Ete3_maker_10.ipynb```](https://github.com/Asklepiad/BI_project_2022/blob/main/scripts/README.md#ete3_maker_10py).
