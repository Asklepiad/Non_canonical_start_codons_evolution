# Non-canonical start-codons and where to find them
Project about evolution of non-canonical start-codons in bacteria and connections between gene function and its start-codon.

## Short description

## Table of content

[Short description](https://github.com/Asklepiad/BI_project_2022/tree/main#short-description)

[Structure of repository](https://github.com/Asklepiad/BI_project_2022/tree/main#structure-of-repository)

[Introduction](https://github.com/Asklepiad/BI_project_2022/tree/main#introduction)

[Aims and objectives](https://github.com/Asklepiad/BI_project_2022/tree/main#aims-and-objectives)

[Data](https://github.com/Asklepiad/BI_project_2022/tree/main#data)

[Workflow](https://github.com/Asklepiad/BI_project_2022/tree/main#workflow)

[Workflow plan and technical properties](https://github.com/Asklepiad/BI_project_2022/tree/main#workflow-plan-and-technical-properties)

[Installation](https://github.com/Asklepiad/BI_project_2022/tree/main#installation)

[Options](https://github.com/Asklepiad/BI_project_2022/tree/main#options)

[Example input and output](https://github.com/Asklepiad/BI_project_2022/tree/main#example-input-and-output)

[Additional info](https://github.com/Asklepiad/BI_project_2022/tree/main#additional-info)

[Results and discussion](https://github.com/Asklepiad/BI_project_2022/tree/main#results-and-discussion)

[Future plans](https://github.com/Asklepiad/BI_project_2022/tree/main#future-plans)

[Literature](https://github.com/Asklepiad/BI_project_2022/tree/main#literature)

[Links](https://github.com/Asklepiad/BI_project_2022/tree/main#links)

[Authors](https://github.com/Asklepiad/BI_project_2022/tree/main#authors)

[Feedback](https://github.com/Asklepiad/BI_project_2022/tree/main#feedback)

## Structure of repository
> Data

All data for reproducing of our investigation and json file for test run are situated in `data` folder. You can read more detailed description of the folder content in `data` folder [readme]()  (*in development*).

> Scripts

Conda environment creator for pipeline, master script, all components of pipeline, some additional scripts and requirements file are located in `scripts` folder. You can read more detailed description of the folder content in `scripts` folder [readme]()  (*in development*).

> Test output

The test running on *Salinibacter ruber* results are situated in `test output` folder (*in development*). You can read more detailed description of the folder content in  `test output` folder [readme]()  (*in development*).

> Illustrations

The images for illustrating repositry's readme files are located in `illustrations` folder (*in development*). You can read more detailed description of the folder content in  `illustrations` folder [readme]()  (*in development*).

## Introduction

It is known procariotes have non only ATG as start-codon. There are about 10% genes in bacterial genome (more strict proportion depends on taxon) with non-canonical start-codons (NCSCs): GTG, TTG and some other in negligibly small proportion. There are not very much data about causes of existance and "gene preferencies" of NCSCs. The GTG and TTG are usually binds weaker with the ribosome relative to canonical ATG start-codone.

## Aims and objectives

1. Find the connection between function of genes and canonicity of its start-codone.
2. Evaluate distribution of bacterial NCSCs from different taxons and ecological niches, find the difference (if it exists) between different groups of bacteria by their function-start-codons patterns.
3. Create pipeline for analysing bacteria.
4. Run pipeline on many organisms from different taxons and ecological niches.

## Data

First stage of investigation was done on thirty-four bacterial species. Complete list of species is available on `data/all_species.txt` file. Species were selected by material of D. Nikolaeva article preprint about generalists and specialists (*article haven't published yet, therefore we can't post any reference or link*).

The json file with precomputed links is situated in `data/id_lists.json`. This file is a main part complete pipeline's input.

## Workflow

The principal scheme of pipeline is on the figure below. ![pipeline of investigation](./illustrations/Pipeline.png)

### Workflow plan and technical properties

#### Technical properties
> The pipeline was run on the local machine with OS Ubuntu 20.04. Package manager was conda 23.3.1 (miniconda). Python version was 3.9.13. R version was 4.2.2 Bash version was 5.0.17(1)-release (x86_64-pc-linux-gnu).

> The pipeline was run on the server with CentOS Linux release 8.5.2111. Package manager was conda 23.3.1 (anaconda). Python version in appropriate virtual environment was 3.9.13. R version was 4.2.2 Bash version was 4.4.20(1)-release (x86_64-redhat-linux-gnu). Workload manager was SLURM 20.11.9-3.1

#### Workflow plan
 There are thirteen scripts in pipeline properly (including two versions of master-script) and some scripts for preparing data.
 
##### Master scripts

There are two master scripts in pipeline for two modes of computing. 
- The more flexible variant is ```folder_creators_server0.sh```. It creates all appropriate directories and subdirectories for one bacteria in input. Use script ```server_executor_start_codons.sbatch``` for the running commands on servers with slurm.
- The more common way is using ```json.sh```.
- In addition you can use parts of pipeline separately. 

##### Pipeline parts properly

1. ```Parsing_NCBI_1.py``` 
2. ```prokka_annotate2.sh``` **(сделать ссылки)**
3. ```First_table_creating3.py```
4. ```proteinortho_script4.sh``` **(сделать ссылки)**
5. ```Muscle_preparing_5.py```
6. ```Statscript.R``` 
7. ```muscle_align6.sh``` **(сделать ссылки)** or prank_align6.sh **(сделать ссылки)**
8. ```modeltree_maker7.sh``` **(сделать ссылки)**
9. ```raxml_tree8.sh``` **(сделать ссылки)**
10. ```Ete3_maker_10.py```

##### Additional scripts

- ```check_strains.py``` finds "bad" assemblies with many genes in ortological rows with gaps on start positions. Script will be implemented in the extended version of pipeline, as part of "believer" module.

- ```download_strains_before.py``` downloads genomes from NCBI, dividing organisms by number of assemblies. The script's work's result is an input for a panakota pipeline **(сделать ссылки)**, which excludes evolutionary close organisms.

- ```download_strains_after.py``` processed output of panakota pipeline, returns json file with dictionary, where key is bacteria's name, and value is a list of links.
- ```environments_script_01.sh``` creates conda virtual environments. If you want to have another virtual environments, or want to set them manually, or run all programms from base conda environment (we highly recommend not to do that), you need to rewrite running scripts and environments names in```folder_creators_server0.sh``` .

### Installation

For the start of work you need to clone this repository on the local maschine or on the server. 
Firstly you need to have conda on your maschine. You can check existance of conda, if you will type `conda` in the command line. If you haven't conda, download it from the [version archive](https://repo.anaconda.com/archive/). You need to choose your OS for downloading. We recommend you to download 


1. You need to clone the repository.
 1. Touch the button `code` on the upper right corner of the screen.
 2. Copy https address.
 3. Open command line with bash.
 4. Change directory, using `cd` command, until you open directory, where you want to copy repository.
 5. Type ```git clone <https you have copied on the previous stage>```
2. After repository had been cloned type ```cd BI_project_2022/scripts``` and after then ```./environments_script_01.sh```
3. When dialog messages will appeared on the screen, you need to push `y` button and enter. It will occurs some times.
4. You are ready for the analysis.

### Options

### Example input and output

### Additional info

### Problems and errors

## Results and discussion

## Future plans

1. Run pipeline on bigger amount of bacteria for verification hypotheses, created on this stage of investigation.
2. Improve quality control and belief in finding position of start-codon (creating "belivier" module with some checkpoint and data filtering points).
3. Found candidat genes for experimental verification in wet lab.

## Literature

1. Belinky, F., Rogozin, I. B., & Koonin, E. V. (2017). Selection on start codons in prokaryotes and potential compensatory nucleotide substitutions. Scientific reports, 7(1), 12422. https://doi.org/10.1038/s41598-017-12619-6
2. Cao, X., & Slavoff, S. A. (2020). Non-AUG start codons: Expanding and regulating the small and alternative ORFeome. Experimental cell research, 391(1), 111973. https://doi.org/10.1016/j.yexcr.2020.111973
3. Galperin, M. Y., Kristensen, D. M., Makarova, K. S., Wolf, Y. I., & Koonin, E. V. (2019). Microbial genome analysis: the COG approach. Briefings in bioinformatics, 20(4), 1063–1070. https://doi.org/10.1093/bib/bbx117
4. Moldovan, M. A., & Gelfand, M. S. (2018). Pangenomic Definition of Prokaryotic Species and the Phylogenetic Structure of Prochlorococcus spp. Frontiers in microbiology, 9, 428. https://doi.org/10.3389/fmicb.2018.00428
5. Costa SS, Guimarães LC, Silva A, Soares SC, Baraúna RA. First Steps in the Analysis of Prokaryotic Pan-Genomes. Bioinformatics and Biology Insights. 2020;14. doi:10.1177/1177932220938064

## Links

## Authors

 - Rubinova Valeria Sergeevna (Bioinformatics institute)
 - Sotnikov Bogdan Vladimirovich (Kyrgyz-Russian Slavic university, Bioinformatics institute)
 - Bochkareva Olga Olegovna* (Institute of Science and Technology Austria)

\* project supervisor 

## Feedback


