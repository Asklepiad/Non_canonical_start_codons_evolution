#!/usr/bin/Rscript --vanilla
print(getwd())

install.packages("ggplot2",repos = "http://cran.us.r-project.org")
install.packages("dplyr",repos = "http://cran.us.r-project.org")
install.packages("fastR2",repos = "http://cran.us.r-project.org")
install.packages("gmodels",repos = "http://cran.us.r-project.org")
install.packages("stringr",repos = "http://cran.us.r-project.org")
install.packages("argparser",repos = "http://cran.us.r-project.org")
install.packages("glue",repos = "http://cran.us.r-project.org")
install.packages("reshape2",repos = "http://cran.us.r-project.org")
install.packages("car",repos = "http://cran.us.r-project.org")
install.packages("hash",repos = "http://cran.us.r-project.org")
#package_installer <- function(package){
#  if (!require(package, character.only=T, quietly=T)) {
#    install.packages(package)
#    library(package, character.only=T)
#  }else{library(package, character.only=T)}
#}
#lapply(c("ggplot2",
#         "dplyr",
#         "fastR2",
#         "gmodels",
#         "stringr",
#         "argparser",
#         "glue",
#         "reshape2",
#         "car"), package_installer)



package_installer <- function(package){
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos = "http://cran.us.r-project.org")
    library(package, character.only=T)
  }else{library(package, character.only=T)}
}
lapply(c("ggplot2",
         "dplyr",
         "fastR2",
         "gmodels",
         "stringr",
         "argparser",
         "glue",
         "reshape2",
         "car",
         "hash"), package_installer)

library(ggplot2)
library(car)
library(gmodels)
library(stringr)
library(argparser)
library(glue)
library(reshape2)
library(dplyr)
library(fastR2)
library(hash)

parser <- arg_parser("Takes folder name")
parser <- add_argument(parser, arg="org_short", type="character", help="Folder name (short organism name)")
arguments = parse_args(parser)
org_short = arguments$org_short

# Uploading the data ####

#org_short <- "B_bifidum"  # For starting from RStudio
print(getwd())
path <- glue("../{org_short}/data/")
setwd(path)
print(getwd())
options(scipen = 999)

summary_rows = read.csv("summary_rows_prokka.csv")
start_codons2 = read.csv("start_codons2_prokka.csv")
summary_rows <- as_tibble(summary_rows)
start_codons2 <- as_tibble(start_codons2)


# Defining se function ####
se <- function(vector){
  n = length(vector)
  return(sd(vector,na.rm=T)/(n**(1/2)))
}


# Common statistics
all_sc_distr <- table(summary_rows$start_codone)


# Computing core-shell-cloude ####
summary_rows$gene_group = "NA"
max_Strain = max(summary_rows$Species)
summary_rows$gene_group = sapply(summary_rows$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
start_codons2$gene_group = sapply(start_codons2$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
# Fixing factor variables ####
summary_rows$gene_group <- as.factor(summary_rows$gene_group)
summary_rows$uniformity <- as.factor(summary_rows$uniformity)
summary_rows$start_codone <- as.factor(summary_rows$start_codone)
summary_rows$type_of_DNA_source <- as.factor(summary_rows$type_of_DNA_source)
summary_rows$p_c_unity <- as.factor(summary_rows$p_c_unity)
#summary_rows$new_sc <- as.factor(summary_rows$new_sc)
#summary_rows$new_length = as.numeric(lapply(summary_rows$new_als, str_length))
# Dividing to three subsets by gene group    ####
core_genes = subset(summary_rows, gene_group=="core")
shell_genes = subset(summary_rows, gene_group=="shell")
cloud_genes = subset(summary_rows, gene_group=="cloud")
cshc_num <- table(summary_rows$gene_group)

# Distributions of sc in cshc ####
abs_core <- table(core_genes$start_codone)
abs_shell <- table(shell_genes$start_codone)
abs_cloud <- table(cloud_genes$start_codone)
core_sc_distr <- round(prop.table(abs_core),3)
shell_sc_distr <- round(prop.table(abs_shell),3)
cloud_sc_distr <- round(prop.table(abs_cloud),3)
# Percent of uniform non-canonical start-codones 
core_nc_unif <- round((sum(subset(start_codons2, gene_group=="core" & uniformity=="same" & ATG==0)$Genes))/sum(subset(start_codons2, gene_group=="core")$Genes),3)
shell_nc_unif <- round((sum(subset(start_codons2, gene_group=="shell" & uniformity=="same" & ATG==0)$Genes))/sum(subset(start_codons2, gene_group=="shell")$Genes),3)
cloud_nc_unif <- round((sum(subset(start_codons2, gene_group=="cloud" & uniformity=="same" & ATG==0)$Genes))/sum(subset(start_codons2, gene_group=="cloud")$Genes),3)
# Percent of ortologus rows with at least one nc start-codone
core_or_wnc <- round((nrow(subset(start_codons2, gene_group=="core" & ATG<Species)))/nrow(subset(start_codons2, gene_group=="core")),3)
shell_or_wnc <- round((nrow(subset(start_codons2, gene_group=="shell" & ATG<Species)))/nrow(subset(start_codons2, gene_group=="shell")),3)
cloud_or_wnc <- round((nrow(subset(start_codons2, gene_group=="cloud" & ATG<Species)))/nrow(subset(start_codons2, gene_group=="cloud")),3)

# Some other basic statistics ####
unif_abs <- length(summary_rows$uniformity[summary_rows$uniformity == "same"])
unif_perc <- round(unif_abs/length(summary_rows$uniformity), 3)
core_abs_unif <- core_genes %>% 
  filter(uniformity=="same") %>% 
  nrow
shell_abs_unif <- shell_genes %>% 
  filter(uniformity=="same") %>% 
  nrow
cloud_abs_unif <- cloud_genes %>% 
  filter(uniformity=="same") %>% 
  nrow
core_perc_unif <- round(core_abs_unif/length(core_genes$uniformity), 3)
shell_perc_unif <- round(shell_abs_unif/length(shell_genes$uniformity), 3)
cloud_perc_unif <- round(cloud_abs_unif/length(cloud_genes$uniformity), 3)
# U-curves ####
## Computing type of start-codone for ortologus row
start_codons2$start_type <- as.factor(ifelse(start_codons2$uniformity == "different", "different", # For colorising the U-curve
                                   ifelse(start_codons2$ATG > 0, "ATG",
                                   ifelse(start_codons2$GTG > 0, "GTG", "TTG"))))
uc_wd <- ggplot(start_codons2, aes(x=Species, fill=start_type))+
  geom_bar()+
  labs(x="Species in ortologus row",
       y="Number of rows")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=24,face="bold"))
ggsave(glue("../figures/{org_short}_uc_wd.png"),  width = 30, height = 20, units = "cm", dpi = 700)
row_starts <- table(start_codons2$start_type)

# boxplots ####
for_bp_of_all <- summary_rows %>% 
  group_by(p_c_unity) %>% 
  summarise(ATG=sum(start_codone=="ATG"),
            GTG=sum(start_codone=="GTG"),
            TTG=sum(start_codone=="TTG")) %>% 
  melt(.) %>% 
  mutate(specie=org_short) %>% 
  rename(start_codon=variable,
         number=value)
write.csv(for_bp_of_all, glue("./{org_short}_for_common_boxplot.csv"))
all_sc_distr_bp <- ggplot()+
  geom_boxplot(data=for_bp_of_all, aes(x=specie, y=number, color=start_codon))
ggsave(glue("../figures/{org_short}_boxplot_all.png"),  width = 30, height = 20, units = "cm", dpi = 700)

violins_scs <- ggplot(summary_rows)+
  geom_violin(aes(x=start_codone, y=Species, fill=start_codone))
ggsave(glue("../figures/{org_short}_violins_scs.png"),  width = 30, height = 20, units = "cm", dpi = 700)
boxplots_scs <- ggplot(summary_rows)+
  geom_boxplot(aes(x=start_codone, y=Species, color=start_codone))
print(getwd())
ggsave(glue("../figures/{org_short}_boxplots_scs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

# Sc distributions per assembly
pc_levels <- levels(summary_rows$p_c_unity)


# common data ####
common_tables <- lapply(pc_levels, function(x) prop.table(table(subset(summary_rows, p_c_unity==x)$start_codone)))
common_atg <- sapply(common_tables, function(x) x["ATG"])
common_gtg <- sapply(common_tables, function(x) x["GTG"])
common_ttg <- sapply(common_tables, function(x) x["TTG"])
###common_other <- sapply(common_tables, function(x) x["Other"])
common_atg_stats <- c(mean(common_atg, na.rm=T), mean(common_atg, na.rm=T)-1.96*se(common_atg), mean(common_atg, na.rm=T)+1.96*se(common_atg))
common_gtg_stats <- c(mean(common_gtg, na.rm=T), mean(common_gtg, na.rm=T)-1.96*se(common_gtg), mean(common_gtg, na.rm=T)+1.96*se(common_gtg))
common_ttg_stats <- c(mean(common_ttg, na.rm=T), mean(common_ttg, na.rm=T)-1.96*se(common_ttg), mean(common_ttg, na.rm=T)+1.96*se(common_ttg))
###common_other_stats <- c(mean(common_other, na.rm=T), mean(common_other, na.rm=T)-1.96*se(common_other), mean(common_other, na.rm=T)+1.96*se(common_other))

# core data
core_tables <- lapply(pc_levels, function(x) prop.table(table(subset(summary_rows, p_c_unity==x & gene_group=="core")$start_codone)))
core_atg <- sapply(core_tables, function(x) x["ATG"])
core_gtg <- sapply(core_tables, function(x) x["GTG"])
core_ttg <- sapply(core_tables, function(x) x["TTG"])
###core_other <- sapply(core_tables, function(x) x["Other"])
core_atg_stats <- c(mean(core_atg, na.rm=T), mean(core_atg, na.rm=T)-1.96*se(core_atg), mean(core_atg, na.rm=T)+1.96*se(core_atg))
core_gtg_stats <- c(mean(core_gtg, na.rm=T), mean(core_gtg, na.rm=T)-1.96*se(core_gtg), mean(core_gtg, na.rm=T)+1.96*se(core_gtg))
core_ttg_stats <- c(mean(core_ttg, na.rm=T), mean(core_ttg, na.rm=T)-1.96*se(core_ttg), mean(core_ttg, na.rm=T)+1.96*se(core_ttg))
###core_other_stats <- c(mean(core_other, na.rm=T), mean(core_other, na.rm=T)-1.96*se(core_other), mean(core_other, na.rm=T)+1.96*se(core_other))

# shell data
shell_tables <- lapply(pc_levels, function(x) prop.table(table(subset(summary_rows, p_c_unity==x & gene_group=="shell")$start_codone)))
shell_atg <- sapply(shell_tables, function(x) x["ATG"])
shell_gtg <- sapply(shell_tables, function(x) x["GTG"])
shell_ttg <- sapply(shell_tables, function(x) x["TTG"])
###shell_other <- sapply(shell_tables, function(x) x["Other"])
shell_atg_stats <- c(mean(shell_atg, na.rm=T), mean(shell_atg, na.rm=T)-1.96*se(shell_atg), mean(shell_atg, na.rm=T)+1.96*se(shell_atg))
shell_gtg_stats <- c(mean(shell_gtg, na.rm=T), mean(shell_gtg, na.rm=T)-1.96*se(shell_gtg), mean(shell_gtg, na.rm=T)+1.96*se(shell_gtg))
shell_ttg_stats <- c(mean(shell_ttg, na.rm=T), mean(shell_ttg, na.rm=T)-1.96*se(shell_ttg), mean(shell_ttg, na.rm=T)+1.96*se(shell_ttg))
###shell_other_stats <- c(mean(shell_other, na.rm=T), mean(shell_other, na.rm=T)-1.96*se(shell_other), mean(shell_other, na.rm=T)+1.96*se(shell_other))

# cloud data
cloud_tables <- lapply(pc_levels, function(x) prop.table(table(subset(summary_rows, p_c_unity==x & gene_group=="cloud")$start_codone)))
cloud_atg <- sapply(cloud_tables, function(x) x["ATG"])
cloud_gtg <- sapply(cloud_tables, function(x) x["GTG"])
cloud_ttg <- sapply(cloud_tables, function(x) x["TTG"])
###cloud_other <- sapply(cloud_tables, function(x) x["Other"])
cloud_atg_stats <- c(mean(cloud_atg, na.rm=T), mean(cloud_atg, na.rm=T)-1.96*se(cloud_atg), mean(cloud_atg, na.rm=T)+1.96*se(cloud_atg))
cloud_gtg_stats <- c(mean(cloud_gtg, na.rm=T), mean(cloud_gtg, na.rm=T)-1.96*se(cloud_gtg), mean(cloud_gtg, na.rm=T)+1.96*se(cloud_gtg))
cloud_ttg_stats <- c(mean(cloud_ttg, na.rm=T), mean(cloud_ttg, na.rm=T)-1.96*se(cloud_ttg), mean(cloud_ttg, na.rm=T)+1.96*se(cloud_ttg))
###cloud_other_stats <- c(mean(cloud_other, na.rm=T), mean(cloud_other, na.rm=T)-1.96*se(cloud_other), mean(cloud_other, na.rm=T)+1.96*se(cloud_other))

# Visualizing cshc vs scs ####
error_bar_df <- as.data.frame(rbind(common_atg_stats, common_gtg_stats, 
                      common_ttg_stats, 
                      core_atg_stats, core_gtg_stats, 
                      core_ttg_stats, 
                      shell_atg_stats, shell_gtg_stats,
                      shell_ttg_stats, 
                      cloud_atg_stats, cloud_gtg_stats, 
                      cloud_ttg_stats))
error_bar_df <- cbind(error_bar_df, rep(c("common", "core", "shell", "cloud"), each=3),
                      rep(c("ATG", "GTG", "TTG")),
                      c("common_atg", "common_gtg", 
                        "common_ttg", 
                        "core_atg", "core_gtg", 
                        "core_ttg",
                        "shell_atg", "shell_gtg",
                        "shell_ttg",
                        "cloud_atg", "cloud_gtg", 
                        "cloud_ttg"))
colnames(error_bar_df) <- c("mean", "lower_ci_bound", "upper_ci_bound", "gene_group", "start_codone", "name")

positions <- c("common_atg", "common_gtg", 
               "common_ttg",
               "core_atg", "core_gtg", 
               "core_ttg",
               "shell_atg", "shell_gtg",
               "shell_ttg",
               "cloud_atg", "cloud_gtg", 
               "cloud_ttg")

# CShC vs SCs plot ####
cshc_scs <- ggplot(error_bar_df)+    
  geom_pointrange(aes(x=name, y=mean, ymin=lower_ci_bound, ymax=upper_ci_bound, group=start_codone, color=start_codone))+
  scale_x_discrete(limits = positions)+
  theme(axis.text.x = element_text(angle = 45))
ggsave(glue("../figures/{org_short}_CShC_scs_eb.png"),  width = 30, height = 20, units = "cm", dpi = 700)

# Uniformity and related information ####
unif_distr <- prop.table(table(summary_rows$uniformity))
unif_nc <- subset(summary_rows, uniformity == "same" & start_codone != "ATG")
unif_nc_table <- table(unif_nc$product)   # List of the genes with uniform non-canonical genes and their frequencies

list_ncs <- unique(as.data.frame(unif_nc_table))
write.csv(list_ncs, glue("./{org_short}_noncanonic_products.csv"))
write.csv(list_ncs, glue("./{org_short}_noncanonic_products.csv"))

# COG statistics ####
# hash
cog_hash <- hash()
cog_abbreviations <- LETTERS
cog_descriptions <- c("rna_proc_and_mod", "chromatin", "energy", "cell_cycle", 
                     "aminoacid", "nucleotide", "carbohydrate", 
                     "coensime", "lipid", "translation_and_ribosomes",  
                     "transcription", "repl_reco_repa",
                     "cell_wall", "cell_motility", "protein_posttrans", "inorganic",  
                     "secondary_metabolites", "general_function_only",
                     "unknown", "signal_transduction", "vesiculs_and_secretion", 
                     "defense", "extracel", "mobilome", "nuclear structure", "cytoskeleton")
sapply(1:26, function(x) cog_hash[[cog_abbreviations[x]]] <- cog_descriptions[x])

## Creating tibble
cog_columns_all <- summary_rows %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1))                      # Choosing all
cog_columns_atg <- summary_rows %>%                          # Choosing atg
  filter(start_codone == "ATG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1))
cog_columns_gtg <- summary_rows %>%                          # Choosing gtg
  filter(start_codone == "GTG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1))
cog_columns_ttg <- summary_rows %>%                          # Choosing ttg
  filter(start_codone == "TTG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1))
# Creating function
cog_names <- sapply(colnames(cog_columns_all), function(x) cog_hash[[x]])
#cog_names <- c("unknown", "transcription", "cell_cycle", "aminoacid", "inorganic", "motility", 
#               "carbohydrate", "lipid", "protein_posttrans", "translation_and_ribosomes", 
#               "mobilome", "cytosceleton", "secondary_metabolites", "vesiculs_and_secretion",
#               "extracel", "chromatine", "general_function_only", "defense", "nucleotide", 
#               "rna_proc_and_mod",  "energy", "cell_wall", "signal_transduction", "coensime", "repl_reco_repa")

#cog_names <- c("translation_and_ribosomes", "rna_proc_and_mod", "transcription", "repl_reco_repa", "chromatine", 
#               "cell_cycle", "defense", "signal_transduction", "cell_wall", "motility", 
#               "cytosceleton", "extracel", "vesiculs_and_secretion", "protein_posttrans", "mobilome",
#               "energy", "carbohydrate", "aminoacid", "nucleotide", "coensime", 
#               "lipid", "inorganic", "secondary_metabolites", "general_function_only", "unknown")
cog_stat_fun <- function(tibble){
  cog_stat <- sapply(tibble, 
         function(x) round(sum(x == 1)/nrow(tibble), 5))
  cog_stat_wilson <- sapply(tibble, 
                            function(x) round(wilson.ci(sum(x == 1), nrow(tibble)),5))
  cog_stat_df <- as.data.frame(cbind(cog_stat, t(cog_stat_wilson)))
  return(cog_stat_df)
}

## Creating pivot table
cog_stat_all <- cog_stat_fun(cog_columns_all)
cog_stat_atg <- cog_stat_fun(cog_columns_atg)
cog_stat_gtg <- cog_stat_fun(cog_columns_gtg)
cog_stat_ttg <- cog_stat_fun(cog_columns_ttg)
cog_pivot <- cbind(cog_names, cog_stat_all, cog_stat_atg, cog_stat_gtg, cog_stat_ttg)   # Summarising table
colnames(cog_pivot) <- c("cog_names", "cog_stat_all", "cog_min_all", "cog_max_all",
                         "cog_stat_atg", "cog_min_atg", "cog_max_atg",
                         "cog_stat_gtg", "cog_min_gtg", "cog_max_gtg",
                         "cog_stat_ttg", "cog_min_ttg", "cog_max_ttg")
write.csv(cog_columns_all, glue("./{org_short}_C_psittaci_COG_precomputed.csv"))

## Drawing errorbars
colors <- c("ATG" = "red", "GTG" = "dark green", "TTG" = "blue", "all" = "black")
ggplot(cog_pivot)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_all, ymin=cog_min_all, ymax=cog_max_all, color="all"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_atg, ymin=cog_min_atg, ymax=cog_max_atg, color="ATG"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_gtg, ymin=cog_min_gtg, ymax=cog_max_gtg, color="GTG"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_ttg, ymin=cog_min_ttg, ymax=cog_max_ttg, color="TTG"))+
  labs("color" = "Legend")+
  scale_color_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 45))

cog_pivot_without_s_and_nulls <- cog_pivot %>% 
  filter(cog_names != "S",
         cog_stat_all != 0)

write.csv(cog_pivot_without_s_and_nulls, file=glue("./{org_short}_C_psittaci_COG_start.csv"))

cog_sc_eb <- ggplot(cog_pivot)+   # Потом сохранить в переменную с именем, как после графика
  geom_pointrange(aes(x=cog_names, y=cog_stat_all, ymin=cog_min_all, ymax=cog_max_all, color="all"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_atg, ymin=cog_min_atg, ymax=cog_max_atg, color="ATG"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_gtg, ymin=cog_min_gtg, ymax=cog_max_gtg, color="GTG"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_ttg, ymin=cog_min_ttg, ymax=cog_max_ttg, color="TTG"), alpha=1)+
  labs("color" = "Legend")+
  xlab(label = "COG categories")+
  ylab(label = "COG frequency of all genes")+
  scale_color_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 35),
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        axis.title.x = element_text(vjust = 13))
  #scale_x_discrete(labels=c("rna_proc_and_mod", "energy", "cell_cycle", 
  #                          "aminoacid", "nucleotide", "carbohydrate", 
   #                         "coensime", "lipid", "translation_and_ribosomes",
    #                        "transcription", "repl_reco_repa", 
     #                       "cell_wall", "inorganic", "protein_posttrans", 
      #                      "motility", "secondary_metabolites", "general_function_only",
       #                     "unknown", "signal_transduction", "vesiculs_and_secretion", 
        #                    "defense", "extracel", "mobilome"))
cog_sc_eb
ggsave(glue("../figures/{org_short}_cog_sc_eb.png"),  width = 30, height = 20, units = "cm", dpi = 700)


cog_sc_eb_short <- ggplot(cog_pivot_without_s_and_nulls)+   # Потом сохранить в переменную с именем, как после графика
  geom_pointrange(aes(x=cog_names, y=cog_stat_all, ymin=cog_min_all, ymax=cog_max_all, color="all"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_atg, ymin=cog_min_atg, ymax=cog_max_atg, color="ATG"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_gtg, ymin=cog_min_gtg, ymax=cog_max_gtg, color="GTG"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_ttg, ymin=cog_min_ttg, ymax=cog_max_ttg, color="TTG"), alpha=1)+
  labs("color" = "Legend")+
  xlab(label = "COG")+
  ylab(label = "COG frequency of all genes")+
  scale_color_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=24,face="bold"))
  #scale_x_discrete(labels=c("rna_proc_and_mod", "energy", "cell_cycle", 
   #                         "aminoacid", "nucleotide", "carbohydrate", 
    #                        "coensime", "lipid", "translation_and_ribosomes",
     #                       "transcription", "repl_reco_repa", 
      #                      "cell_wall", "inorganic", "protein_posttrans", 
       #                     "motility", "secondary_metabolites", "general_function_only",
        #                    "signal_transduction", "vesiculs_and_secretion", "defense", 
         #                   "extracel", "mobilome"))
cog_sc_eb_short
ggsave(glue("../figures/{org_short}_cog_sc_eb_short.png"),  width = 30, height = 20, units = "cm", dpi = 700)
### COG formal test (exact Fisher) #####
#### Function for COG formal ######
cog_formal_atg <- sapply(cog_columns_atg, 
                     function(x) sum(x == 1))
cog_formal_gtg <- sapply(cog_columns_gtg, 
                         function(x) sum(x == 1))
cog_formal_ttg <- sapply(cog_columns_ttg, 
                         function(x) sum(x == 1))
cog_formal <- as.data.frame(cbind(cog_formal_atg,
                                  cog_formal_gtg,
                                  cog_formal_ttg))
## Excluded rows with 0 values in a row
cog_formal_without_zeros <- cog_formal %>%     # Initialisation of shorted dataset
  filter(cog_formal_atg+cog_formal_gtg+cog_formal_ttg > 0)
write.csv(cog_formal_without_zeros, glue("{org_short}_cogs.csv"))
#fi <- fisher.test(cog_formal_without_zeros, simulate.p.value = TRUE, B=150000)   # Computings for absolute values
#hi <- chisq.test(cog_formal_without_zeros)
#hi$expected
#mosaicplot(cog_formal_without_zeros, color=T, shade=T, xlab="Start-codon", ylab="Function")
#cfwz <- mutate(rowwise(cog_formal_without_zeros), total = sum(c_across(1:3)))    
#
#cfwz_percents <- cfwz %>%   # Dataset with realtive data (percents)
#  transmute(atg = cog_formal_atg/total,
#            gtg = cog_formal_gtg/total,
#            ttg = cog_formal_ttg/total)
#chisq.test(cfwz_percents)   # One bad idea

#b <- cog_formal_without_zeros %>%
#  table() %>%

#d <- summary_rows %>% 
#  select(start_codone, J) %>% 
#  table %>% 
#  CrossTable(prop.r = FALSE, prop.c = FALSE, prop.t = FALSE, prop.chisq = FALSE, fisher = TRUE, simulate.p.value=TRUE)


cfwz_percents <- cog_formal_without_zeros %>%   # Dataset with realtive data (percents)
  transmute(atg = (cog_formal_atg/sum(cog_formal_atg)),
            gtg = (cog_formal_gtg/sum(cog_formal_gtg)),
            ttg = (cog_formal_ttg/sum(cog_formal_ttg)),
            total = (cog_formal_atg + cog_formal_gtg + cog_formal_ttg)/(sum(cog_formal_atg) + sum(cog_formal_gtg + cog_formal_ttg)))
#chisq.test(cfwz_percents)   # One bad idea
#cfwz_percents <- t(cfwz_percents)
#sapply(c(1:3), function(x) chisq.test(rbind(cfwz_percents[x, ], cfwz_percents[4, ])))


## Part of non-canonic-starts in orto-rows
atg_content <- ggplot(start_codons2)+
  geom_point(aes(x=Species, y=ATG), alpha=0.1, color="red")
ggsave(glue("../figures/{org_short}_atg_content.png"),  width = 30, height = 20, units = "cm", dpi = 700)
gtg_content <- ggplot(start_codons2)+
  geom_point(aes(x=Species, y=GTG), alpha=0.1, color="green")
ggsave(glue("../figures/{org_short}_gtg_content.png"),  width = 30, height = 20, units = "cm", dpi = 700)
ttg_content <- ggplot(start_codons2)+
  geom_point(aes(x=Species, y=TTG), alpha=0.1, color="blue")
ggsave(glue("../figures/{org_short}_ttg_content.png"),  width = 30, height = 20, units = "cm", dpi = 700)
## Cogs in non-canons per assembly
assembly_list <- unique(summary_rows$p_c_unity)
cog_list <- colnames(cog_columns_all)

smth <- nrow(subset(summary_rows, p_c_unity==0 & S==1 & start_codone!="ATG"))/nrow(subset(summary_rows, p_c_unity==0 & S==1))


# Length of non-uniform ####
## Are there any difference between uniform and non-uniform length diff into ortorologus rows
# normalised_median_deviation <- function(vector){
#   med = median(vector)
#   square_distances = sapply(vector, function(x) (med - x)**2)
#   result = ((sum(square_distances)/(length(vector - 1)))**(1/2))/med
#   return(result)
# }
# 
# # Choosing non-uniform genes
# non_unif_summary_rows <- summary_rows %>% 
#   filter(uniformity == "different")
# 
# # Counting mdn
# non_uni <- non_unif_summary_rows %>% 
#   group_by(ortologus_row) %>% 
#   summarise(max_length=max(length),
#             med_length=median(length),
#             mdn=normalised_median_deviation(length))
# View(non_uni)
# 
# # Choosing uniform genes
# unif_summary_rows <- summary_rows %>% 
#   filter(uniformity == "same", Genes>1)
# 
# # Counting mdn
# 
# uni <- unif_summary_rows %>% 
#   group_by(ortologus_row) %>% 
#   summarise(max_length=max(length),
#             med_length=median(length),
#             mdn=normalised_median_deviation(length))
# View(uni)
# 
# shapiro.test(non_uni$mdn)
# shapiro.test(uni$mdn)
# # The data is strongly unnormal distributed
# median(non_uni$mdn)
# median(uni$mdn)
# # How many percent of rows have deviation more than 10%
# length(non_uni$mdn[non_uni$mdn>0.1])/length(non_uni$mdn)
# length(non_uni$mdn[uni$mdn>0.1])/length(uni$mdn)
# # How many percent of rows haven't deviation
# length(non_uni$mdn[non_uni$mdn==0])/length(non_uni$mdn)
# length(non_uni$mdn[uni$mdn==0])/length(uni$mdn)
# # Formal test
# wilcox.test(non_uni$mdn, uni$mdn)
# 
# ## Is removing of short sequences removed diversity too?
# pre_non <- non_unif_summary_rows %>% 
#   group_by(ortologus_row) %>% 
#   summarise(max_codon=(max(table(start_codone))/sum(table(start_codone))),
#             codone_name_pre=names(table(start_codone))[which.max(table(start_codone))])
# post_non <- non_unif_summary_rows %>%
#   group_by(ortologus_row) %>%
#   filter(length>=median(length)) %>% 
#   summarise(max_codon=(max(table(start_codone))/sum(table(start_codone))),
#             codone_name_post=names(table(start_codone))[which.max(table(start_codone))])
# diversity_changing <- left_join(pre_non, post_non, by="ortologus_row") %>% 
#   mutate(difference=max_codon.y-max_codon.x) %>% 
#   mutate(codone_differs=ifelse(codone_name_pre==codone_name_post, "No", "Yes"))
# 
# 
# # Alternative start-codons
# 
# ## Changing in length after new alignment
# # Choosing non-uniform genes
# non_unif_summary_rows <- summary_rows %>% 
#   filter(uniformity == "different")
# 
# # Counting mdn
# non_uni <- non_unif_summary_rows %>% 
#   group_by(ortologus_row) %>% 
#   summarise(max_length_before=max(length),
#             med_length_before=median(length),
#             mdn_before=normalised_median_deviation(length),
#             max_length_after=max(new_length),
#             med_length_after=median(new_length),
#             mdn_after=normalised_median_deviation(new_length),
#             mdn_change=mdn_before-mdn_after)
# View(non_uni)
# 
# # Comparing maximal frequent start-codons per row before and after length increasing, and after correcting
# pre_post_corr_non <- non_unif_summary_rows %>% 
#   group_by(ortologus_row) %>% 
#   summarise(max_codon_pre=(max(table(start_codone))/sum(table(start_codone))),
#             codone_name_pre=names(table(start_codone))[which.max(table(start_codone))],
#             max_codon_post=(max(table(new_sc))/sum(table(new_sc))),
#             codone_name_post=names(table(new_sc))[which.max(table(new_sc))],
#             max_codon_corr=(max(table(new_sc_c))/sum(table(new_sc_c))),
#             codone_name_corr=names(table(new_sc_c))[which.max(table(new_sc_c))],
#             change_pre_post=max_codon_post-max_codon_pre,
#             change_pre_corr=max_codon_corr-max_codon_pre,
#             is_diff_pre_post=(codone_name_pre!=codone_name_post),
#             is_diff_pre_corr=(codone_name_pre!=codone_name_corr))
# View(pre_post_corr_non)
# 
# ggplot(summary_rows)+
#   geom_histogram(aes(x=delta), binwidth = 50)
# 
# Gene group proportion visualisation
prop_gene_group <- start_codons2 %>% 
  group_by(gene_group, start_type) %>%
  summarise(count = n()) %>%
  filter(gene_group != "NA" & is.na(start_type)==FALSE)

# Barplots CShC - number ####
positions = c("core", "shell", "cloud")
or_bar_abs <-  ggplot(prop_gene_group) +
  scale_x_discrete(limits = positions)+
  geom_col(aes(x = gene_group, y = count, fill = start_type))
ggsave(glue("../figures/{org_short}_or_bar_abs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

or_bar_rel <- ggplot(prop_gene_group, aes(x = gene_group, y = count, fill = start_type)) +
  geom_bar(stat="identity", position="fill")+
  scale_x_discrete(limits = positions)+
  xlab(label = "Pangenome fraction")+
  ylab(label = "Percent of ortologus rows with same start-codons")+
  annotate("text", x=c(1, 2, 3, 1, 2, 0.9, 2, 3, 1.1, 3), 
           y = c(0.75, 0.75, 0.75, 0.17, 0.45, 0.06, 0.17, 0.17, 0.01, 0.05), 
           label=c("612", "9", "293", "201", "5", "43", "7", "54", "17", "34"),
           size = 12)+
  theme(axis.text=element_text(size=20, face="bold"),
        axis.title.x=element_text(size=24,face="bold"),
        axis.title.y=element_text(size=24,face="bold"),
        legend.title = element_text(size=18),
        legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=14))+
  scale_fill_discrete(name = "Start codon")
ggsave(glue("../figures/{org_short}_or_bar_rel.png"),  width = 30, height = 20, units = "cm", dpi = 700)


half_blood_gene <- summary_rows %>%
  filter(ortologus_row==297)
#half_blood_gene$product


# Lera ####  
# Сделать нормальные подписи и прочие украшательства
distr_scs_common = as.data.frame(table(summary_rows$start_codone))
genes_col_abs <- ggplot(distr_scs_common, aes(x=Var1, y=Freq, fill=Var1))+
  geom_bar(stat="identity")
ggsave(glue("../figures/{org_short}_genes_col_abs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

scs_perrow <- start_codons2 %>% 
  select(ATG:TTG) %>% 
  apply(., 1, function(x) ifelse(x[1] >= x[2] & x[1] >= x[3], "ATG",
                                 ifelse(x[2] >= x[3] & x[2] > x[1], "GTG", "TTG")))
or_col_abs <- ggplot(as.data.frame(table(scs_perrow)), aes(x=scs_perrow, y=Freq, fill=scs_perrow))+
  geom_bar(stat="identity")
ggsave(glue("../figures/{org_short}_OR_col_abs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

start_codons2$scs_perrow <- scs_perrow
U_curve_wod <- ggplot(start_codons2, aes(x=Species, fill=scs_perrow))+
  geom_bar()+
  labs(x="Species in ortologus row",
       y="Number of rows")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=24,face="bold"))
ggsave(glue("../figures/{org_short}_UC_wod.png"),  width = 30, height = 20, units = "cm", dpi = 700)

cog_abs_nc <- summary_rows %>%    # Сколько и каких когов (абс) у неканоник
  filter(start_codone != "ATG") %>% 
  select(S:L) %>% 
  summarise_all(sum)

cog_perc_nc <- summary_rows %>%    # Сколько и каких когов (%) у неканоник
  filter(start_codone != "ATG") %>% 
  select(S:L) %>% 
  summarise_all(mean)

cog_abs_atg <- summary_rows %>%         # Сколько и каких когов (абс) у каноник
  filter(start_codone == "ATG") %>% 
  select(S:L) %>% 
  summarise_all(sum)   

cog_perc_atg <- summary_rows %>%         # Сколько и каких когов (%) у каноник
  filter(start_codone == "ATG") %>% 
  select(S:L) %>% 
  summarise_all(mean)  

cog_abs_gtg <- summary_rows %>%         # Сколько и каких когов (абс) у GTG
  filter(start_codone == "GTG") %>% 
  select(S:L) %>% 
  summarise_all(sum)   

cog_perc_gtg <- summary_rows %>%         # Сколько и каких когов (%) у GTG
  filter(start_codone == "GTG") %>% 
  select(S:L) %>% 
  summarise_all(mean)  

cog_abs_ttg <- summary_rows %>%         # Сколько и каких когов (абс) у TTG
  filter(start_codone == "TTG") %>% 
  select(S:L) %>% 
  summarise_all(sum)   

cog_perc_ttg <- summary_rows %>%         # Сколько и каких когов (%) у TTG
  filter(start_codone == "TTG") %>% 
  select(S:L) %>% 
  summarise_all(mean)  

cog_abs_all <- summary_rows %>%         # Сколько и каких когов (абс) у всех
  select(S:L) %>% 
  summarise_all(sum)   

cog_perc_all <- summary_rows %>%         # Сколько и каких когов (%) у всех
  select(S:L) %>% 
  summarise_all(mean)

cog_abs_perc <- as.data.frame(rbind(cog_abs_nc, cog_perc_nc, cog_abs_atg, cog_perc_atg, 
      cog_abs_gtg, cog_perc_gtg, cog_abs_ttg, 
      cog_perc_ttg, cog_abs_all, cog_perc_all))
rownames(cog_abs_perc) <- c("cog_abs_nc", "cog_perc_nc", "cog_abs_atg", "cog_perc_atg", 
                            "cog_abs_gtg", "cog_perc_gtg", "cog_abs_ttg", 
                            "cog_perc_ttg", "cog_abs_all", "cog_perc_all")
cog_abs_perc <- t(cog_abs_perc)
write.csv(cog_abs_perc, glue("./{org_short}_cog_stat_per_sc.csv"))

have_cogs_absolute <- summary_rows %>%         # Сколько генов с когом S
  select(S) %>% 
  table()
names(have_cogs_absolute) <- c("Have cogs", "Haven't cogs")

rows_wcogs <- summary_rows %>%    # Число рядов, которым не приписалось ни одного кога
  group_by(ortologus_row) %>% 
  summarise("cognot" = sum(S)/n()) %>% 
  filter(cognot == 1) %>% 
  count()

rows_numbers <- start_codons2 %>% 
  filter(Species >= Genes) %>% 
  nrow

# Lavrenty idea ####

ass_cog_atg <- summary_rows %>% 
  filter(start_codone == "ATG") %>% 
  select(p_c_unity, length:ortologus_row) %>% 
  select(p_c_unity, 3:(ncol(.)-1)) %>% 
  group_by(p_c_unity) %>% 
  summarize_all(sum)

ass_cog_gtg <- summary_rows %>% 
  filter(start_codone == "GTG") %>% 
  select(p_c_unity, length:ortologus_row) %>% 
  select(p_c_unity, 3:(ncol(.)-1)) %>%
  group_by(p_c_unity) %>% 
  summarize_all(sum)

ass_cog_ttg <- summary_rows %>% 
  filter(start_codone == "TTG") %>% 
  select(p_c_unity, length:ortologus_row) %>% 
  select(p_c_unity, 3:(ncol(.)-1)) %>%
  group_by(p_c_unity) %>% 
  summarize_all(sum)

cogs_names <- summary_rows %>%
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>%
  summarise_all(sum) %>% 
  t %>% 
  as.data.frame %>% 
  filter(V1>0) %>% 
  rownames

coeffs <- c(table(summary_rows$start_codone)[1]/table(summary_rows$start_codone)[3], table(summary_rows$start_codone)[1]/table(summary_rows$start_codone)[2])
ass_sc_creator <- function(cog){
  ass_sc <- as.data.frame(cbind(pull(ass_cog_atg[cog]), pull(ass_cog_gtg[cog]), pull(ass_cog_ttg[cog])))
  colnames(ass_sc) <- c("ATG", "GTG", "TTG")
  ass_sc$GTG <- ass_sc$GTG*coeffs[2]
  ass_sc$TTG <- ass_sc$TTG*coeffs[1]
  return(ass_sc)
}
coeffs <- c(table(summary_rows$start_codone)[1]/table(summary_rows$start_codone)[3], table(summary_rows$start_codone)[2]/table(summary_rows$start_codone)[3])
ass_sc_creator <- function(cog){
  ass_sc <- as.data.frame(cbind(pull(ass_cog_atg[cog]), pull(ass_cog_gtg[cog]), pull(ass_cog_ttg[cog])))
  colnames(ass_sc) <- c("ATG", "GTG", "TTG")
  ass_sc$GTG <- round(ass_sc$GTG/coeffs[2],5)
  ass_sc$ATG <- round(ass_sc$ATG/coeffs[1],5)
  return(ass_sc)
}
ass_sc_meltor <- function(ass_sc){
  ass_sc <- melt(ass_sc)
  names(ass_sc) <- c("start_codone", "count")
  return(ass_sc)
}
cog_posthoc <- function(cog_df){
  AG <- wilcox.test(cog_df$ATG, cog_df$GTG, exact=FALSE)
  AT <- wilcox.test(cog_df$ATG, cog_df$TTG, exact=FALSE)
  GT <- wilcox.test(cog_df$TTG, cog_df$GTG, exact=FALSE)
  us <- c(AG$statistic[[1]], AT$statistic[[1]], GT$statistic[[1]])
  ps <- c(AG$p.value, AT$p.value, GT$p.value)
  aps <- p.adjust(ps, method = "bonferroni")
  res_df <- data.frame(us, ps, aps)
  rownames(res_df) <- c("ATG_GTG", "ATG_TTG", "GTG_TTG")
  return(res_df)
}

ass_sc_cog_pre_list <- lapply(cogs_names, ass_sc_creator)
names(ass_sc_cog_pre_list) <- cogs_names
ass_sc_cog_list <- lapply(ass_sc_cog_pre_list, ass_sc_meltor)
kw_cogs_data <- lapply(ass_sc_cog_list, function(x) kruskal.test(count ~ start_codone, data = x))
kw_cogs_pval <- sapply(kw_cogs_data, function(x) x$p.value)
bh_corr <- p.adjust(kw_cogs_pval, method = "bonferroni")
bh_corr_v <- bh_corr[bh_corr<0.05]
mwu_b_list <- lapply(ass_sc_cog_pre_list, cog_posthoc)
mwu_df <- as.data.frame(do.call("rbind", mwu_b_list))
write.csv(mwu_df, glue("./{org_short}_cog_sc_mwu.csv"))


# Report ####  
report_part_basic <- glue("Basic statistics:
                     Number of genes = {nrow(summary_rows)}
                     Number of ortologus rows = {nrow(start_codons2)}
                     Number of ortologus rows without paralogs = {rows_numbers}
                     All genes start-codon distribution: ATG {all_sc_distr[1]}
                                                         GTG {all_sc_distr[2]}
                                                         TTG {all_sc_distr[3]}")

report_cshc_1 <- glue("Core, shell and cloud (partitions):
                     There are {cshc_num[2]} genes in core.
                     Core genes' start-codon distribution: ATG {core_sc_distr[1]} 
                                                           GTG {core_sc_distr[2]} 
                                                           TTG {core_sc_distr[3]}
                     There are {cshc_num[4]} genes in shell.                                      
                     Shell genes' start-codon distribution: ATG {shell_sc_distr[1]} 
                                                            GTG {shell_sc_distr[2]} 
                                                            TTG {shell_sc_distr[3]}
                     There are {cshc_num[1]} genes in cloud.
                     Cloud genes' start-codon distribution: ATG {cloud_sc_distr[1]} 
                                                            GTG {cloud_sc_distr[2]} 
                                                            TTG {cloud_sc_distr[3]}")


report_cshc_2 <- glue("Core, shell and cloud (absolute):
                     Core genes' start-codon number: ATG {abs_core[1]} 
                                                     GTG {abs_core[2]} 
                                                     TTG {abs_core[3]}
                     Shell genes' start-codon number: ATG {abs_shell[1]} 
                                                      GTG {abs_shell[2]} 
                                                      TTG {abs_shell[3]}
                     Cloud genes' start-codon number: ATG {abs_cloud[1]} 
                                                      GTG {abs_cloud[2]} 
                                                      TTG {abs_cloud[3]}")

report_unif <- glue("Unformity of ortologus rows:
                    There are {unif_abs} uniform rows. It is {unif_perc*100} percents of all ortologus rows.
                    Core genes contains {core_abs_unif} uniform rows. It is {core_perc_unif*100} percents.
                    Shell genes contains {shell_abs_unif} uniform rows. It is {shell_perc_unif*100} percents.
                    Cloud genes contains {cloud_abs_unif} uniform rows. It is {cloud_perc_unif*100} percents.
                    Percent of uniform rows with non-canonical start-codon in core is {core_nc_unif*100}
                    Percent of uniform rows with non-canonical start-codon in shell is {shell_nc_unif*100}
                    Percent of uniform rows with non-canonical start-codon in cloud is {cloud_nc_unif*100}")

report_scs_add <- glue("Additional data about start codons:
                       There are {core_or_wnc} percents of core ortologus rows with at least one start-codon
                       There are {shell_or_wnc} percents of shell ortologus rows with at least one start-codon
                       There are {cloud_or_wnc} percents of cloud ortologus rows with at least one start-codon")

report_cogs2 <- glue("The number of rows without cogs (or with cog 'S') is equal {rows_wcogs[[1]]} or {round(rows_wcogs[[1]]/rows_numbers, 3)*100} percents")

write.table(c(report_part_basic, report_cshc_1, report_cshc_2, report_unif, report_scs_add, report_cogs2), glue("{org_short}_statistical_short_report.txt"))


