#!/usr/bin/Rscript --vanilla
print(getwd())

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
         "hash",
         "readr"), package_installer)


# Uploading the data ####

# For starting from RStudio add your bacteria's short name -- first letter of genus, underscore, first nine (or lesser, if hasn't) letters of specie's name.
org_short <- "E_coli"  
print(getwd())
path <- glue("../results/{org_short}")
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
## Computing type of start-codone for ortologuos row
start_codons2$start_type <- as.factor(ifelse(start_codons2$uniformity == "different", "different", # For colorising the U-curve
                                   ifelse(start_codons2$ATG > 0, "ATG",
                                   ifelse(start_codons2$GTG > 0, "GTG", "TTG"))))
uc_wd <- ggplot(start_codons2, aes(x=Species, fill=start_type))+
  geom_bar()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=24,face="bold"))+
  scale_fill_manual(values = c("ATG" = "#F8766D",
                               "GTG" = "#00BA38",
                               "TTG" = "#619CFF",
                               "different" = "#CD9600"))+
  labs("fill" = "Start codon of ortologus row",
       x="Species in ortologuos row",
       y="Number of rows")
ggsave(glue("./figures/{org_short}_uc_wd.png"),  width = 30, height = 20, units = "cm", dpi = 700)
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
write.csv(for_bp_of_all, glue("./data/{org_short}_for_common_boxplot.csv"))
all_sc_distr_bp <- ggplot()+
  geom_boxplot(data=for_bp_of_all, aes(x=specie, y=number, color=start_codon))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))+
  labs("color" = "Start codon",
       x = "Specie",
       y = "Number of genes")
ggsave(glue("./figures/{org_short}_boxplot_all.png"),  width = 30, height = 20, units = "cm", dpi = 700)

violins_scs <- ggplot(summary_rows)+
  geom_violin(aes(x=start_codone, y=Species, fill=start_codone))+
  theme(axis.text=element_text(size=18),
      axis.title=element_text(size=24,face="bold"))+
  labs("color" = "Start codon",
       x = "Start codon",
       y = "Strains")
ggsave(glue("./figures/{org_short}_violins_scs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

boxplots_scs <- ggplot(summary_rows)+
  geom_boxplot(aes(x=start_codone, y=Species, color=start_codone))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))+
  labs("color" = "Start codon",
       x = "Start codon",
       y = "Strains")
print(getwd())
ggsave(glue("./figures/{org_short}_boxplots_scs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

# Sc distributions per assembly
pc_levels <- levels(summary_rows$p_c_unity)


# common data ####
common_tables <- lapply(pc_levels, function(x) prop.table(table(subset(summary_rows, p_c_unity==x)$start_codone)))
common_atg <- sapply(common_tables, function(x) x["ATG"])
common_gtg <- sapply(common_tables, function(x) x["GTG"])
common_ttg <- sapply(common_tables, function(x) x["TTG"])

common_atg_stats <- c(mean(common_atg, na.rm=T), mean(common_atg, na.rm=T)-1.96*se(common_atg), mean(common_atg, na.rm=T)+1.96*se(common_atg))
common_gtg_stats <- c(mean(common_gtg, na.rm=T), mean(common_gtg, na.rm=T)-1.96*se(common_gtg), mean(common_gtg, na.rm=T)+1.96*se(common_gtg))
common_ttg_stats <- c(mean(common_ttg, na.rm=T), mean(common_ttg, na.rm=T)-1.96*se(common_ttg), mean(common_ttg, na.rm=T)+1.96*se(common_ttg))

# core data
core_tables <- lapply(pc_levels, function(x) prop.table(table(subset(summary_rows, p_c_unity==x & gene_group=="core")$start_codone)))
core_atg <- sapply(core_tables, function(x) x["ATG"])
core_gtg <- sapply(core_tables, function(x) x["GTG"])
core_ttg <- sapply(core_tables, function(x) x["TTG"])

core_atg_stats <- c(mean(core_atg, na.rm=T), mean(core_atg, na.rm=T)-1.96*se(core_atg), mean(core_atg, na.rm=T)+1.96*se(core_atg))
core_gtg_stats <- c(mean(core_gtg, na.rm=T), mean(core_gtg, na.rm=T)-1.96*se(core_gtg), mean(core_gtg, na.rm=T)+1.96*se(core_gtg))
core_ttg_stats <- c(mean(core_ttg, na.rm=T), mean(core_ttg, na.rm=T)-1.96*se(core_ttg), mean(core_ttg, na.rm=T)+1.96*se(core_ttg))

# shell data
shell_tables <- lapply(pc_levels, function(x) prop.table(table(subset(summary_rows, p_c_unity==x & gene_group=="shell")$start_codone)))
shell_atg <- sapply(shell_tables, function(x) x["ATG"])
shell_gtg <- sapply(shell_tables, function(x) x["GTG"])
shell_ttg <- sapply(shell_tables, function(x) x["TTG"])

shell_atg_stats <- c(mean(shell_atg, na.rm=T), mean(shell_atg, na.rm=T)-1.96*se(shell_atg), mean(shell_atg, na.rm=T)+1.96*se(shell_atg))
shell_gtg_stats <- c(mean(shell_gtg, na.rm=T), mean(shell_gtg, na.rm=T)-1.96*se(shell_gtg), mean(shell_gtg, na.rm=T)+1.96*se(shell_gtg))
shell_ttg_stats <- c(mean(shell_ttg, na.rm=T), mean(shell_ttg, na.rm=T)-1.96*se(shell_ttg), mean(shell_ttg, na.rm=T)+1.96*se(shell_ttg))

# cloud data
cloud_tables <- lapply(pc_levels, function(x) prop.table(table(subset(summary_rows, p_c_unity==x & gene_group=="cloud")$start_codone)))
cloud_atg <- sapply(cloud_tables, function(x) x["ATG"])
cloud_gtg <- sapply(cloud_tables, function(x) x["GTG"])
cloud_ttg <- sapply(cloud_tables, function(x) x["TTG"])

cloud_atg_stats <- c(mean(cloud_atg, na.rm=T), mean(cloud_atg, na.rm=T)-1.96*se(cloud_atg), mean(cloud_atg, na.rm=T)+1.96*se(cloud_atg))
cloud_gtg_stats <- c(mean(cloud_gtg, na.rm=T), mean(cloud_gtg, na.rm=T)-1.96*se(cloud_gtg), mean(cloud_gtg, na.rm=T)+1.96*se(cloud_gtg))
cloud_ttg_stats <- c(mean(cloud_ttg, na.rm=T), mean(cloud_ttg, na.rm=T)-1.96*se(cloud_ttg), mean(cloud_ttg, na.rm=T)+1.96*se(cloud_ttg))

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
  theme(axis.text.x = element_text(angle = 45),
    axis.text=element_text(size=18),
    axis.title=element_text(size=24,face="bold"))+
  labs("color" = "Start codon",
       x = "Pangenome group and start codon",
       y = "Mean proportion of start codon type")
ggsave(glue("./figures/{org_short}_CShC_scs_eb.png"),  width = 30, height = 20, units = "cm", dpi = 700)

# Uniformity and related information ####
unif_distr <- prop.table(table(summary_rows$uniformity))
unif_nc <- subset(summary_rows, uniformity == "same" & start_codone != "ATG")
unif_nc_table <- table(unif_nc$product)   # List of the genes with uniform non-canonical genes and their frequencies

list_ncs <- unique(as.data.frame(unif_nc_table))
write.csv(list_ncs, glue("./data/{org_short}_noncanonic_products.csv"))

# COG statistics ####
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
write.csv(cog_columns_all, glue("./data/{org_short}_COG_precomputed.csv"))

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
  filter(cog_names != "unknown",
         cog_stat_all != 0)

write.csv(cog_pivot_without_s_and_nulls, file=glue("./data/{org_short}_COG_start.csv"))

cog_sc_eb <- ggplot(cog_pivot)+   # Потом сохранить в переменную с именем, как после графика
  geom_pointrange(aes(x=cog_names, y=cog_stat_all, ymin=cog_min_all, ymax=cog_max_all, color="all"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_atg, ymin=cog_min_atg, ymax=cog_max_atg, color="ATG"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_gtg, ymin=cog_min_gtg, ymax=cog_max_gtg, color="GTG"), alpha=1)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_ttg, ymin=cog_min_ttg, ymax=cog_max_ttg, color="TTG"), alpha=1)+
  labs("color" = "Start codon")+
  xlab(label = "COG categories")+
  ylab(label = "COG frequency of all genes")+
  scale_color_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 35),
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        axis.title.x = element_text(vjust = 13))
cog_sc_eb
ggsave(glue("./figures/{org_short}_cog_sc_eb.png"),  width = 30, height = 20, units = "cm", dpi = 700)


cog_sc_eb_short <- ggplot(cog_pivot_without_s_and_nulls)+   # Потом сохранить в переменную с именем, как после графика
  geom_pointrange(aes(x=cog_names, y=cog_stat_all, ymin=cog_min_all, ymax=cog_max_all, color="all"), size=1.3, alpha=0.6)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_atg, ymin=cog_min_atg, ymax=cog_max_atg, color="ATG"), size=1.3, alpha=0.6)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_gtg, ymin=cog_min_gtg, ymax=cog_max_gtg, color="GTG"), size=1.3, alpha=0.6)+
  geom_pointrange(aes(x=cog_names, y=cog_stat_ttg, ymin=cog_min_ttg, ymax=cog_max_ttg, color="TTG"), size=1.3, alpha=0.6)+
  labs("color" = "Legend",
       x = "COG",
       y = "Given start-codon's type proportion\nin given COG")+
  scale_color_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 27))+
  theme(axis.text=element_text(size=23),
        axis.title=element_text(size=24,face="bold"),
        axis.title.x = element_text(vjust=12.5),
        legend.position = "none",
        plot.margin = margin(,4,-2,, "cm"))
cog_sc_eb_short
ggsave(glue("./figures/{org_short}_cog_sc_eb_short.png"),  width = 30, height = 20, units = "cm", dpi = 700)
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

cfwz_percents <- cog_formal_without_zeros %>%   # Dataset with realtive data (percents)
  transmute(atg = (cog_formal_atg/sum(cog_formal_atg)),
            gtg = (cog_formal_gtg/sum(cog_formal_gtg)),
            ttg = (cog_formal_ttg/sum(cog_formal_ttg)),
            total = (cog_formal_atg + cog_formal_gtg + cog_formal_ttg)/(sum(cog_formal_atg) + sum(cog_formal_gtg + cog_formal_ttg)))

## Part of non-canonic-starts in orto-rows
atg_content <- ggplot(start_codons2)+
  geom_point(aes(x=Species, y=ATG), alpha=0.1, color="red")+
  theme(axis.text=element_text(size=18),
    axis.title=element_text(size=24,face="bold"))+
  labs(x = "Strains",
       y = "Percent of ATG in ortologus row")
ggsave(glue("./figures/{org_short}_atg_content.png"),  width = 30, height = 20, units = "cm", dpi = 700)

gtg_content <- ggplot(start_codons2)+
  geom_point(aes(x=Species, y=GTG), alpha=0.1, color="green")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))+
  labs(x = "Strains",
       y = "Percent of GTG in ortologus row")
ggsave(glue("./figures/{org_short}_gtg_content.png"),  width = 30, height = 20, units = "cm", dpi = 700)

ttg_content <- ggplot(start_codons2)+
  geom_point(aes(x=Species, y=TTG), alpha=0.1, color="blue")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))+
  labs(x = "Strains",
       y = "Percent of TTG in ortologus row")
ggsave(glue("./figures/{org_short}_ttg_content.png"),  width = 30, height = 20, units = "cm", dpi = 700)

## Cogs in non-canons per assembly
assembly_list <- unique(summary_rows$p_c_unity)
cog_list <- colnames(cog_columns_all)

smth <- nrow(subset(summary_rows, p_c_unity==0 & S==1 & start_codone!="ATG"))/nrow(subset(summary_rows, p_c_unity==0 & S==1))

# Gene group proportion visualisation
prop_gene_group <- start_codons2 %>% 
  group_by(gene_group, start_type) %>%
  summarise(count = n()) %>%
  filter(gene_group != "NA" & is.na(start_type)==FALSE)


# Fisher on steroids ####

# Actual approach
## Sum of all cogs
al_sum <- summary_rows %>% 
  select(length:ortologus_row) %>%
  select(3:ncol(.)-1) %>%
  sum

## Proportion of COGs
prop_cog <- summary_rows %>% 
  select(length:ortologus_row) %>%
  select(3:ncol(.)-1) %>% 
  summarise(across(everything(), function(x) round(sum(x)/al_sum, 4))) %>% 
  select(where(function(x) sum(x)>0))

## Proportion of start codons and creating reference tibble
pre_fisher_ref <- tibble(start_codone = c("ATG", "GTG", "TTG"),
                         count = table(summary_rows$start_codone)) 


## Dataset fo further analysis: x -- COGs, y -- start-codons
pre_fisher <- summary_rows %>% 
  select(start_codone, length:ortologus_row) %>%
  select(1, 4:ncol(.)-1)  %>%  # Why 4? Why not 3?
  group_by(start_codone) %>%
  summarize(across(1:ncol(.)-1, sum))  %>% # Why we need to use "-1"?
  select(start_codone, where(~ is.factor(.x) || (is.numeric(.x) && sum(.x)>0)))



list_sc_num <- list(c(1,2), c(1,3), c(2,3))
vector_cog <- colnames(prop_cog)

fisher_list <- lapply(list_sc_num, 
       function(x) lapply(vector_cog, 
                          function(y) fisher.test(
                            rbind(as.integer((pre_fisher_ref[x,2] * as.double(prop_cog[y]))[1:2,2]),
                                  as.vector(pre_fisher[x, y])[[1]]))))

names(fisher_list) <- c("ATG_GTG", "ATG_TTG", "GTG_TTG")
names(fisher_list[["ATG_GTG"]]) <- vector_cog
names(fisher_list[["ATG_TTG"]]) <- vector_cog
names(fisher_list[["GTG_TTG"]]) <- vector_cog

bonf_corr <- length(vector_cog) * length(list_sc_num)
lapply(c("ATG_GTG", "ATG_TTG", "GTG_TTG"),
       function(x) lapply(vector_cog, 
                          function(y) fisher_list[[x]][[y]]["corrected_p_value"] = <подставить нужное>)

fisher.test(matrix(c(15332, 1228, 16875, 1264), nrow=2))

  

# Gene selecting for eggnog-mapper ####
pre_egg_or <- summary_rows %>% 
  filter(gene_group == "core") %>% 
  select(id, source, aa_sequence, ortologus_row, length)

first_part_egg <- pre_egg_or %>% 
  group_by(ortologus_row) %>% 
  filter(length==median(length)) %>% 
  slice_sample(n=1) %>% 
  ungroup

second_part_egg <- pre_egg_or %>% 
  group_by(ortologus_row) %>% 
  filter(length!=median(length)) %>% 
  distinct(ortologus_row, aa_sequence, .keep_all = TRUE) %>% 
  ungroup

egg_pre_fasta <- rbind(first_part_egg, second_part_egg)

## Fasta creating
egg_fasta <- egg_pre_fasta %>% 
  select(id, source, aa_sequence) %>% 
  mutate(name = paste0(">", paste(id, source, sep="__")))

fasta_prewriter <- function(x, col1, col2){
  cat(x[col1])
  cat("\n")
  cat(x[col2])
  cat("\n")
}

sink(glue("../../data/for_eggnog_{org_short}.fa"))
apply(egg_fasta, 1, function(x) fasta_prewriter(x, "name", "aa_sequence"))
sink()


# Barplots CShC - number ####
positions = c("core", "shell", "cloud")
or_bar_abs <-  ggplot(prop_gene_group) +
  scale_x_discrete(limits = positions)+
  geom_col(aes(x = gene_group, y = count, fill = start_type))+
  scale_fill_manual(values = c("ATG" = "#F8766D",
                               "GTG" = "#00BA38",
                               "TTG" = "#619CFF",
                               "different" = "#CD9600"))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))+
  labs("fill" = "Start codon",
       x = "Gene group",
       y = "Number of ortologus rows")
ggsave(glue("./figures/{org_short}_or_bar_abs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

or_bar_rel <- ggplot(prop_gene_group, aes(x = gene_group, y = count, fill = start_type)) +
  geom_bar(stat="identity", position="fill")+
  scale_x_discrete(limits = positions)+
  xlab(label = "Pangenome fraction")+
  ylab(label = "Percent of ortologus rows \nwith same start-codons")+
  theme(axis.text=element_text(size=20, face="bold"),
        axis.title.x=element_text(size=24,face="bold"),
        axis.title.y=element_text(size=24,face="bold"),
        legend.title = element_text(size=18),
        legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=14))+
  scale_fill_manual(values = c("ATG" = "#F8766D",
                               "GTG" = "#00BA38",
                               "TTG" = "#619CFF",
                               "different" = "#CD9600"))+
  labs("fill" = "Start codon")
ggsave(glue("./figures/{org_short}_or_bar_rel.png"),  width = 30, height = 20, units = "cm", dpi = 700)


half_blood_gene <- summary_rows %>%
  filter(ortologus_row==297)

# Additional statistics ####  

distr_scs_common = as.data.frame(table(summary_rows$start_codone))
genes_col_abs <- ggplot(distr_scs_common, aes(x=Var1, y=Freq, fill=Var1))+
  geom_bar(stat="identity")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))+
  xlab(label = "Start codon")+
  ylab(label = "Number of genes")+
  labs("fill" = "Start codon")
ggsave(glue("./figures/{org_short}_genes_col_abs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

scs_perrow <- start_codons2 %>% 
  select(ATG:TTG) %>% 
  apply(., 1, function(x) ifelse(x[1] >= x[2] & x[1] >= x[3], "ATG",
                                 ifelse(x[2] >= x[3] & x[2] > x[1], "GTG", "TTG")))

or_col_abs <- ggplot(as.data.frame(table(scs_perrow)), aes(x=scs_perrow, y=Freq, fill=scs_perrow))+
  geom_bar(stat="identity")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))+
  xlab(label = "Start codon")+
  ylab(label = "Number of genes")+
  labs("fill" = "Start codon")
ggsave(glue("./figures/{org_short}_OR_col_abs.png"),  width = 30, height = 20, units = "cm", dpi = 700)

start_codons2$scs_perrow <- scs_perrow
U_curve_wod <- ggplot(start_codons2, aes(x=Species, fill=scs_perrow))+
  geom_bar()+
  labs(x="Species in ortologus row",
       y="Number of rows",
       "fill"="Start codon")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=24,face="bold"))+
  scale_fill_manual(values = c("ATG" = "#F8766D",
                               "GTG" = "#00BA38",
                               "TTG" = "#619CFF"))

ggsave(glue("./figures/{org_short}_UC_wod.png"),  width = 30, height = 20, units = "cm", dpi = 700)




cog_abs_nc <- summary_rows %>%    # COG number and types (abs) in non-canonical SCs
  filter(start_codone != "ATG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(sum)

cog_perc_nc <- summary_rows %>%    # COG number and types (%) in non-canonical SCs
  filter(start_codone != "ATG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(mean)

cog_abs_atg <- summary_rows %>%         # COG number and types (abs) in canonical SCs
  filter(start_codone == "ATG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(sum)   

cog_perc_atg <- summary_rows %>%         # COG number and types (%) in canonical SCs
  filter(start_codone == "ATG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(mean)  

cog_abs_gtg <- summary_rows %>%         # COG number and types (abs) in GTG
  filter(start_codone == "GTG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(sum)   

cog_perc_gtg <- summary_rows %>%         # COG number and types (%) in GTG
  filter(start_codone == "GTG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(mean)  

cog_abs_ttg <- summary_rows %>%         # COG number and types (abs) in TTG
  filter(start_codone == "TTG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(sum)   

cog_perc_ttg <- summary_rows %>%         # COG number and types (%) in TTG
  filter(start_codone == "TTG") %>% 
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(mean)  

cog_abs_all <- summary_rows %>%         # COG number and types (abs) in all SCs
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(sum)   

cog_perc_all <- summary_rows %>%         # COG number and types (%) in all SCs
  select(length:ortologus_row) %>% 
  select(2:(ncol(.)-1)) %>% 
  summarise_all(mean)

cog_abs_perc <- as.data.frame(rbind(cog_abs_nc, cog_perc_nc, cog_abs_atg, cog_perc_atg, 
      cog_abs_gtg, cog_perc_gtg, cog_abs_ttg, 
      cog_perc_ttg, cog_abs_all, cog_perc_all))
rownames(cog_abs_perc) <- c("cog_abs_nc", "cog_perc_nc", "cog_abs_atg", "cog_perc_atg", 
                            "cog_abs_gtg", "cog_perc_gtg", "cog_abs_ttg", 
                            "cog_perc_ttg", "cog_abs_all", "cog_perc_all")
cog_abs_perc <- t(cog_abs_perc)
write.csv(cog_abs_perc, glue("./data/{org_short}_cog_stat_per_sc.csv"))

have_cogs_absolute <- summary_rows %>%         # Number of genes withoue known function
  select(S) %>% 
  table()
names(have_cogs_absolute) <- c("Have cogs", "Haven't cogs")

rows_wcogs <- summary_rows %>%    # Number of ortologous rows without cogs
  group_by(ortologus_row) %>% 
  summarise("cognot" = sum(S)/n()) %>% 
  filter(cognot == 1) %>% 
  count()

rows_numbers <- start_codons2 %>% 
  filter(Species >= Genes) %>% 
  nrow

# Formal test for evaluating COG-SC interactions ####

# ass_cog_atg <- summary_rows %>% 
#   filter(start_codone == "ATG") %>% 
#   select(p_c_unity, length:ortologus_row) %>% 
#   select(p_c_unity, 3:(ncol(.)-1)) %>% 
#   group_by(p_c_unity) %>% 
#   summarize_all(sum)
# 
# ass_cog_gtg <- summary_rows %>% 
#   filter(start_codone == "GTG") %>% 
#   select(p_c_unity, length:ortologus_row) %>% 
#   select(p_c_unity, 3:(ncol(.)-1)) %>%
#   group_by(p_c_unity) %>% 
#   summarize_all(sum)
# 
# ass_cog_ttg <- summary_rows %>% 
#   filter(start_codone == "TTG") %>% 
#   select(p_c_unity, length:ortologus_row) %>% 
#   select(p_c_unity, 3:(ncol(.)-1)) %>%
#   group_by(p_c_unity) %>% 
#   summarize_all(sum)
# 
# cogs_names <- summary_rows %>%
#   select(length:ortologus_row) %>% 
#   select(2:(ncol(.)-1)) %>%
#   summarise_all(sum) %>% 
#   t %>% 
#   as.data.frame %>% 
#   filter(V1>0) %>% 
#   rownames
# 
# coeffs <- c(table(summary_rows$start_codone)[1]/table(summary_rows$start_codone)[3], table(summary_rows$start_codone)[1]/table(summary_rows$start_codone)[2])
# ass_sc_creator <- function(cog){
#   ass_sc <- as.data.frame(cbind(pull(ass_cog_atg[cog]), pull(ass_cog_gtg[cog]), pull(ass_cog_ttg[cog])))
#   colnames(ass_sc) <- c("ATG", "GTG", "TTG")
#   ass_sc$GTG <- ass_sc$GTG*coeffs[2]
#   ass_sc$TTG <- ass_sc$TTG*coeffs[1]
#   return(ass_sc)
# }
# coeffs <- c(table(summary_rows$start_codone)[1]/table(summary_rows$start_codone)[3], table(summary_rows$start_codone)[2]/table(summary_rows$start_codone)[3])
# ass_sc_creator <- function(cog){
#   ass_sc <- as.data.frame(cbind(pull(ass_cog_atg[cog]), pull(ass_cog_gtg[cog]), pull(ass_cog_ttg[cog])))
#   colnames(ass_sc) <- c("ATG", "GTG", "TTG")
#   ass_sc$GTG <- round(ass_sc$GTG/coeffs[2],5)
#   ass_sc$ATG <- round(ass_sc$ATG/coeffs[1],5)
#   return(ass_sc)
# }
# ass_sc_meltor <- function(ass_sc){
#   ass_sc <- melt(ass_sc)
#   names(ass_sc) <- c("start_codone", "count")
#   return(ass_sc)
# }
# cog_posthoc <- function(cog_df){
#   AG <- wilcox.test(cog_df$ATG, cog_df$GTG, exact=FALSE)
#   AT <- wilcox.test(cog_df$ATG, cog_df$TTG, exact=FALSE)
#   GT <- wilcox.test(cog_df$TTG, cog_df$GTG, exact=FALSE)
#   us <- c(AG$statistic[[1]], AT$statistic[[1]], GT$statistic[[1]])
#   ps <- c(AG$p.value, AT$p.value, GT$p.value)
#   aps <- p.adjust(ps, method = "bonferroni")
#   res_df <- data.frame(us, ps, aps)
#   rownames(res_df) <- c("ATG_GTG", "ATG_TTG", "GTG_TTG")
#   return(res_df)
# }
# 
# ass_sc_cog_pre_list <- lapply(cogs_names, ass_sc_creator)
# names(ass_sc_cog_pre_list) <- cogs_names
# ass_sc_cog_list <- lapply(ass_sc_cog_pre_list, ass_sc_meltor)
# kw_cogs_data <- lapply(ass_sc_cog_list, function(x) kruskal.test(count ~ start_codone, data = x))
# kw_cogs_pval <- sapply(kw_cogs_data, function(x) x$p.value)
# bh_corr <- p.adjust(kw_cogs_pval, method = "bonferroni")
# bh_corr_v <- bh_corr[bh_corr<0.05]
# mwu_b_list <- lapply(ass_sc_cog_pre_list, cog_posthoc)
# mwu_df <- as.data.frame(do.call("rbind", mwu_b_list))
# write.csv(mwu_df, glue("./data/{org_short}_cog_sc_mwu.csv"))


# Genes of interest ####

# Choosing core orto-rows from organism with ATG proportion lesser then 0.5
summary_rows_ec = read.csv("./results/E_coli/summary_rows_prokka.csv")
start_codons2_ec = read.csv("./results/E_coli/start_codons2_prokka.csv")
max_Strain = max(summary_rows_ec$Species)
summary_rows_ec$gene_group = sapply(summary_rows_ec$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
start_codons2_ec$gene_group = sapply(start_codons2_ec$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
E_coli_interest_rows <- summary_rows_ec %>%
  filter(gene_group == "core") %>%
  group_by(ortologus_row) %>%
  summarize(atg_content = (ATG / Genes)) %>%
  unique %>%
  filter(atg_content < 0.5) %>%
  pull(ortologus_row)
E_coli_interest <- subset(summary_rows_ec, ortologus_row %in% E_coli_interest_rows)
E_coli_interest_product1 <- unique(E_coli_interest$product)
E_coli_interest_product2 <- E_coli_interest %>%
  group_by(product) %>%
  summarize(atg_content = (ATG / Genes)) %>%
  unique %>%
  filter(atg_content < 0.5) %>%
  pull(product)
E_coli_interest_products <- intersect(E_coli_interest_product1, E_coli_interest_product2)
E_coli_interest_products <- E_coli_interest_products[E_coli_interest_products != "hypothetical protein"]

summary_rows_pa = read.csv("./results/P_aeruginos/summary_rows_prokka.csv")
start_codons2_pa = read.csv("./results/P_aeruginos/start_codons2_prokka.csv")
max_Strain = max(summary_rows_pa$Species)
summary_rows_pa$gene_group = sapply(summary_rows_pa$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
start_codons2_pa$gene_group = sapply(start_codons2_pa$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
P_aeruginos_interest_rows <- summary_rows_pa %>%
  filter(gene_group == "core") %>%
  group_by(ortologus_row) %>%
  summarize(atg_content = (ATG / Genes)) %>%
  unique %>%
  filter(atg_content < 0.5) %>%
  pull(ortologus_row)
P_aeruginos_interest <- subset(summary_rows_pa, ortologus_row %in% P_aeruginos_interest_rows)
P_aeruginos_interest_product1 <- unique(P_aeruginos_interest$product)
P_aeruginos_interest_product2 <- P_aeruginos_interest %>%
  group_by(product) %>%
  summarize(atg_content = (ATG / Genes)) %>%
  unique %>%
  filter(atg_content < 0.5) %>%
  pull(product)
P_aeruginos_interest_products <- intersect(P_aeruginos_interest_product1, P_aeruginos_interest_product2)
P_aeruginos_interest_products <- P_aeruginos_interest_products[P_aeruginos_interest_products != "hypothetical protein"]

summary_rows_bc = read.csv("./results/B_cereus/summary_rows_prokka.csv")
start_codons2_bc = read.csv("./results/B_cereus/start_codons2_prokka.csv")
max_Strain = max(summary_rows_bc$Species)
summary_rows_bc$gene_group = sapply(summary_rows_bc$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
start_codons2_bc$gene_group = sapply(start_codons2_bc$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
B_cereus_interest_rows <- summary_rows_bc %>%
  filter(gene_group == "core") %>%
  group_by(ortologus_row) %>%
  summarize(atg_content = (ATG / Genes)) %>%
  unique %>%
  filter(atg_content < 0.5) %>%
  pull(ortologus_row)
B_cereus_interest <- subset(summary_rows_bc, ortologus_row %in% B_cereus_interest_rows)
B_cereus_interest_product1 <- unique(B_cereus_interest$product)
B_cereus_interest_product2 <- B_cereus_interest %>%
  group_by(product) %>%
  summarize(atg_content = (ATG / Genes)) %>%
  unique %>%
  filter(atg_content < 0.5) %>%
  pull(product)
B_cereus_interest_products <- intersect(B_cereus_interest_product1, B_cereus_interest_product2)
B_cereus_interest_products <- B_cereus_interest_products[B_cereus_interest_products != "hypothetical protein"]

intersect(intersect(B_cereus_interest_products, P_aeruginos_interest_products), E_coli_interest_products)
#subset(summary_rows, starts_with())

# Alignment of "ATP synthase subunit b(eta)" ####
b_subunit <- summary_rows_ec |> 
  slice(grep("ATP synthase subunit b", product)) |> 
  select(product, aa_sequence)

b_subunit$product <-  paste(">", gsub(" ", "_", b_subunit$product), 0:99, sep="_")

fasta_prewriter <- function(x, col1, col2){
  cat(x[col1])
  cat("\n")
  cat(x[col2])
  cat("\n")
}

sink("./data/atp_b.fa")
apply(b_subunit, 1, function(x) fasta_prewriter(x, "product", "aa_sequence"))
sink()


       
       
#summary_rows = read.csv("../results/P_aeruginos/summary_rows_prokka.csv")
#start_codons2 = read.csv("../results/P_aeruginos/start_codons2_prokka.csv")
#max_Strain = max(summary_rows$Species)
#start_codons2$gene_group = sapply(start_codons2$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
#P_aeruginos_interest_rows <- unique(subset(start_codons2, gene_group == "core" & ATG/(ATG + GTG + TTG) < 0.33)$cog)
#P_aeruginos_interest_cogs <- P_aeruginos_interest_cogs[P_aeruginos_interest_cogs != "absent"]
#P_aeruginos_interest_products <- unique(subset(summary_rows, cog %in% P_aeruginos_interest_cogs)$product)

#summary_rows = read.csv("../results/B_cereus/summary_rows_prokka.csv")
#start_codons2 = read.csv("../results/B_cereus/start_codons2_prokka.csv")
#max_Strain = max(summary_rows$Species)
#start_codons2$gene_group = sapply(start_codons2$Species, function(x) ifelse(x==max_Strain, "core", ifelse(x==1, "cloud", ifelse((round(max_Strain*0.7,0)>=x) & (x>=round(max_Strain*0.3,0)), "shell", "NA"))))
#B_cereus_interest_rows <- unique(subset(start_codons2, gene_group == "core" & ATG/(ATG + GTG + TTG) < 0.33)$cog)
#B_cereus_interest_cogs <- B_cereus_interest_cogs[B_cereus_interest_cogs != "absent"]
#B_cereus_interest_products <- unique(subset(summary_rows, cog %in% B_cereus_interest_cogs)$product)

consensus_interest_cogs1 <- intersect(intersect(P_aeruginos_interest_products, B_cereus_interest_products), E_coli_interest_products)

summary_rows = read.csv("../results/E_coli/summary_rows_prokka.csv")
lapply(consensus_interest_cogs1, function(x) prop.table(table(subset(summary_rows, product == x)$start_codone)))

summary_rows = read.csv("../results/P_aeruginos/summary_rows_prokka.csv")
lapply(consensus_interest_cogs1, function(x) prop.table(table(subset(summary_rows, product == x)$start_codone)))

summary_rows = read.csv("../results/B_cereus/summary_rows_prokka.csv")
lapply(consensus_interest_cogs1, function(x) prop.table(table(subset(summary_rows, product == x)$start_codone)))

summary_rows %>%
  filter(product )

#Generation the order of COG names
#cog_list <- summary_rows %>%  
#  select(length:ortologus_row) %>% 
#  select(2:(ncol(.)-1)) %>%
#  colnames(.)

# Generation of list with start-codon "weight" per COG
#cog_start_prop_list <- lapply(cog_list, function(x) prop.table(table(subset(summary_rows, summary_rows[x]==1)$start_codone)))
#names(cog_start_prop_list) <- cog_list

#orto_cores <- subset(start_codons2, gene_group == "core")

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

setwd("../../scripts/")


