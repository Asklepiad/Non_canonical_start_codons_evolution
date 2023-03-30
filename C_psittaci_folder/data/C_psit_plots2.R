# Uploading the data ####
library(ggplot2)
library(dplyr)
library(fastR2)
library(stringr)
setwd("~/bioinf/start_codons/BI_project_2022/C_psittaci_folder/C_psittaci/data/")
options(scipen = 999)
#c_psit2 <- read.csv("./c_psittaci2.proteinortho.tsv", sep = "\t")
#View(c_psit2)
summary_rows = read.csv("summary_rows_new.csv")
start_codons2 = read.csv("start_codons2_prokka.csv")
summary_rows <- as_tibble(summary_rows)
start_codons2 <- as_tibble(start_codons2)


# Defining se function ####
se <- function(vector){
  n = length(vector)
  return(sd(vector,na.rm=T)/(n**(1/2)))
}


# Computing core-shell-cloude ####
summary_rows$gene_group = "NA"
max_strain = max(summary_rows$Strains)
summary_rows$gene_group = sapply(summary_rows$Strains, function(x) ifelse(x==max_strain, "core", ifelse(x==1, "cloud", ifelse((round(max_strain*0.7,0)>=x) & (x>=round(max_strain*0.3,0)), "shell", "NA"))))
start_codons2$gene_group = sapply(start_codons2$Strains, function(x) ifelse(x==max_strain, "core", ifelse(x==1, "cloud", ifelse((round(max_strain*0.7,0)>=x) & (x>=round(max_strain*0.3,0)), "shell", "NA"))))
# Fixing factor variables ####
summary_rows$gene_group <- as.factor(summary_rows$gene_group)
summary_rows$type_of_the_gene <- as.factor(summary_rows$type_of_the_gene)
summary_rows$uniformity <- as.factor(summary_rows$uniformity)
summary_rows$start_codone <- as.factor(summary_rows$start_codone)
summary_rows$type_of_DNA_source <- as.factor(summary_rows$type_of_DNA_source)
summary_rows$p_c_unity <- as.factor(summary_rows$p_c_unity)
summary_rows$new_sc <- as.factor(summary_rows$new_sc)
summary_rows$new_length = as.numeric(lapply(summary_rows$new_als, str_length))
# Dividing to three subsets by gene group
core_genes = subset(summary_rows, gene_group=="core")
shell_genes = subset(summary_rows, gene_group=="shell")
cloud_genes = subset(summary_rows, gene_group=="cloud")


# Distributions of sc in different groups ####
core_sc_distr <- prop.table(table(core_genes$start_codone))
shell_sc_distr <- prop.table(table(shell_genes$start_codone))
cloud_sc_distr <- prop.table(table(cloud_genes$start_codone))
# Percent of uniform non-canonical start-codones 
core_nc_unif <- (sum(subset(start_codons2, gene_group=="core" & uniformity=="uniform" & ATG==0)$Genes))/sum(subset(start_codons2, gene_group=="core")$Genes)
shell_nc_unif <- (sum(subset(start_codons2, gene_group=="shell" & uniformity=="uniform" & ATG==0)$Genes))/sum(subset(start_codons2, gene_group=="shell")$Genes)
cloud_nc_unif <- (sum(subset(start_codons2, gene_group=="cloud" & uniformity=="uniform" & ATG==0)$Genes))/sum(subset(start_codons2, gene_group=="cloud")$Genes)
# Percent of ortologus rows with at least one nc start-codone
core_or_wnc <- (nrow(subset(start_codons2, gene_group=="core" & ATG<Strains)))/nrow(subset(start_codons2, gene_group=="core"))
shell_or_wnc <- (nrow(subset(start_codons2, gene_group=="shell" & ATG<Strains)))/nrow(subset(start_codons2, gene_group=="shell"))
cloud_or_wnc <- (nrow(subset(start_codons2, gene_group=="cloud" & ATG<Strains)))/nrow(subset(start_codons2, gene_group=="cloud"))


# U-curves ####
## Computing type of start-codone for ortologus row
start_codons2$start_type <- as.factor(ifelse(start_codons2$uniformity == "non-uniform", "different", # For colorising the U-curve
                                   ifelse(start_codons2$ATG > 0, "ATG",
                                   ifelse(start_codons2$GTG > 0, "GTG", "TTG"))))
ggplot(start_codons2, aes(x=Strains, fill=start_type))+
  geom_bar()+
  labs(x="Strains in ortologus row",
       y="Number of rows")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=24,face="bold"))
table(start_codons2$start_type)

# boxplots ####
ggplot(summary_rows)+
  geom_boxplot(aes(x=start_codone, y=Strains, color=start_codone))
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

# Visualizing
error_bar_df <- rbind(common_atg_stats, common_gtg_stats, 
                      common_ttg_stats, 
                      core_atg_stats, core_gtg_stats, 
                      core_ttg_stats, 
                      shell_atg_stats, shell_gtg_stats,
                      shell_ttg_stats, 
                      cloud_atg_stats, cloud_gtg_stats, 
                      cloud_ttg_stats)
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
error_bar_df <- as.data.frame(error_bar_df)
positions <- c("common_atg", "common_gtg", 
               "common_ttg",
               "core_atg", "core_gtg", 
               "core_ttg",
               "shell_atg", "shell_gtg",
               "shell_ttg",
               "cloud_atg", "cloud_gtg", 
               "cloud_ttg")
error_bar_df <- as.data.frame(error_bar_df)
ggplot(error_bar_df)+
  geom_pointrange(aes(x=name, y=mean, ymin=lower_ci_bound, ymax=upper_ci_bound, group=start_codone, color=start_codone))+
  scale_x_discrete(limits = positions)+
  theme(axis.text.x = element_text(angle = 45),
        #axis.text.y=element_blank())
  )
# Uniformity and related information ####
unif_distr <- prop.table(table(summary_rows$uniformity))
unif_nc <- subset(summary_rows, uniformity == "uniform" & start_codone != "ATG")
unif_nc_table <- table(unif_nc$product)   # List of the genes with uniform non-canonical genes and their frequencies


# COG statistics ####
## Creating tibble
cog_columns_all <- select(summary_rows, S:L)                 # Choosing all
cog_columns_atg <- summary_rows %>%                          # Choosing atg
  filter(start_codone == "ATG") %>% 
  select(S:L)
cog_columns_gtg <- summary_rows %>%                          # Choosing gtg
  filter(start_codone == "GTG") %>% 
  select(S:L)
cog_columns_ttg <- summary_rows %>%                          # Choosing ttg
  filter(start_codone == "TTG") %>% 
  select(S:L)
# Creating function
cog_names <- c("unknown", "transcription", "cell_cycle", "aminoacid", "inorganic", "motility", 
               "carbohydrate", "lipid", "protein_posttrans", "translation_and_ribosomes", 
               "mobilome", "cytosceleton", "secondary_metabolites", "vesiculs_and_secretion",
               "extracel", "chromatine", "general_function_only", "defense", "nucleotide", 
               "rna_proc_and_mod",  "energy", "cell_wall", "signal_transduction", "coensime", "repl_reco_repa")

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
write.csv(cog_columns_all, "./C_psittaci_COG_precomputed.csv")

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

write.csv(cog_pivot_without_s_and_nulls, file="./C_psittaci_COG_start.csv")

present_cog_plot <- ggplot(cog_pivot_without_s_and_nulls)+   # Потом сохранить в переменную с именем, как после графика
  geom_pointrange(aes(x=cog_names, y=cog_stat_all, ymin=cog_min_all, ymax=cog_max_all, color="all"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_atg, ymin=cog_min_atg, ymax=cog_max_atg, color="ATG"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_gtg, ymin=cog_min_gtg, ymax=cog_max_gtg, color="GTG"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_ttg, ymin=cog_min_ttg, ymax=cog_max_ttg, color="TTG"))+
  labs("color" = "Legend")+
  xlab(label = "COG categories")+
  ylab(label = "COG frequency of all genes")+
  scale_color_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 35),
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        axis.title.x = element_text(vjust = 13))+
  annotate("text", x = c(2, 2, 3,  15, 16, 19, 19, 20, 20), 
           y = c(0.07, 0.013, 0.035, 0.055, 0.04, 0.05, 0.01, 0.2, 0.07), 
           label = c("GTG carbohydrate", "TTG carbohydrate", 
                     "GTG cell cycle", 
                     "GTG protein posttranslational modifications",
                     "ATG replication, recombination, reparation",
                     "GTG transcription", "TTG transcription",
                     "TTG translation", "GTG translation"),
           colour = c("dark green", "blue", "dark green",
                      "dark green", "red", "dark green", "blue", 
                      "blue", "dark green"),
           size=7)+
  annotate("rect", xmin = c(1.5, 2.5, 14.5, 15.5, 18.5, 19.5), 
           xmax = c(2.5, 3.5, 15.5, 16.5, 19.5, 20.5), ymin = 0, ymax = 0.21,
             alpha = .2, fill="yellow", colour="magenta")
present_cog_plot
ggsave("./cog_categories.png",  width = 30, height = 15, units = "cm")




present_cog_plot_prev <- ggplot(cog_pivot_without_s_and_nulls)+   # Потом сохранить в переменную с именем, как после графика
  geom_pointrange(aes(x=cog_names, y=cog_stat_all, ymin=cog_min_all, ymax=cog_max_all, color="all"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_atg, ymin=cog_min_atg, ymax=cog_max_atg, color="ATG"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_gtg, ymin=cog_min_gtg, ymax=cog_max_gtg, color="GTG"))+
  geom_pointrange(aes(x=cog_names, y=cog_stat_ttg, ymin=cog_min_ttg, ymax=cog_max_ttg, color="TTG"))+
  labs("color" = "Legend")+
  xlab(label = "COG")+
  ylab(label = "COG frequency of all genes")+
  scale_color_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 45))+
  #annotate("text", x = c(2, 2, 3, 4, 5, 6.9, 8, 8.1, 10, 11, 14, 15, 16, 19, 19, 20, 20), 
   #        y = c(0.07, 0.013, 0.035, 0.015, 0.022, 0.025, 0.053, 0.019, 0.023, 0.052, 0.007, 0.055, 0.04, 0.05, 0.01, 0.2, 0.07), 
    #       label = c("GTG carbohydrate", "TTG carbohydrate", 
     #                "GTG cell cycle", "TTG cell wall", 
      #               "GTG coensime", "GTG defense",
       #              "TTG energy", "GTG energy",
        #             "GTG general functions", "TTG inorganic",
         #            "GTG nucleotide",
          #           "GTG protein posttranslational modifications",
           #          "ATG replication, recombination, reparation",
            #         "GTG transcription", "TTG transcription",
             #        "TTG translation", "GTG translation"),
           #colour = c("darkgreen", "blue", "green", "blue",
            #          "green", "green", "blue", "green",
             #         "green", "blue", "green", "green", 
              #        "red", "green", "blue", "blue", "green"),
          # size=7)+
  #annotate("rect", xmin = c(1.5, 2.5, 14.5, 15.5, 18.5, 19.5), 
   #        xmax = c(2.5, 3.5, 15.5, 16.5, 19.5, 20.5), ymin = 0, ymax = 0.21,
    #       alpha = .2, fill="yellow", colour="magenta")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=24,face="bold"))
present_cog_plot_prev
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
cog_formal_without_zeros <- cog_formal %>% 
  filter(cog_formal_atg+cog_formal_gtg+cog_formal_ttg > 0)
hihi <- chisq.test(cog_formal_without_zeros)
fisher.test(cog_formal_without_zeros)


## Part of non-canonic-starts in orto-rows
ggplot(start_codons2)+
  geom_point(aes(x=Strains, y=ATG), alpha=0.1, color="red")
ggplot(start_codons2)+
  geom_point(aes(x=Strains, y=GTG), alpha=0.1, color="green")
ggplot(start_codons2)+
  geom_point(aes(x=Strains, y=TTG), alpha=0.1, color="blue")

## Cogs in non-canons per assembly
assembly_list <- unique(summary_rows$p_c_unity)
cog_list <- colnames(cog_columns_all)
sapply(assembly_list, function(x) sapply(cog_list, function(y) nrow(subset(summary_rows, p_c_unity==x & y==1 & start_codone!="ATG"))/nrow(subset(summary_rows, p_c_unity==x & y==1))))


nrow(subset(summary_rows, p_c_unity==0 & S==1 & start_codone!="ATG"))/nrow(subset(summary_rows, p_c_unity==0 & S==1))


# Length of non-uniform ####
## Are there any difference between uniform and non-uniform length diff into ortorologus rows
normalised_median_deviation <- function(vector){
  med = median(vector)
  square_distances = sapply(vector, function(x) (med - x)**2)
  result = ((sum(square_distances)/(length(vector - 1)))**(1/2))/med
  return(result)
}

# Choosing non-uniform genes
non_unif_summary_rows <- summary_rows %>% 
  filter(uniformity == "non-uniform")

# Counting mdn
non_uni <- non_unif_summary_rows %>% 
  group_by(ortologus_row) %>% 
  summarise(max_length=max(length),
            med_length=median(length),
            mdn=normalised_median_deviation(length))
View(non_uni)

# Choosing uniform genes
unif_summary_rows <- summary_rows %>% 
  filter(uniformity == "uniform", Genes>1)

# Counting mdn

uni <- unif_summary_rows %>% 
  group_by(ortologus_row) %>% 
  summarise(max_length=max(length),
            med_length=median(length),
            mdn=normalised_median_deviation(length))
View(uni)

shapiro.test(non_uni$mdn)
shapiro.test(uni$mdn)
# The data is strongly unnormal distributed
median(non_uni$mdn)
median(uni$mdn)
# How many percent of rows have deviation more than 10%
length(non_uni$mdn[non_uni$mdn>0.1])/length(non_uni$mdn)
length(non_uni$mdn[uni$mdn>0.1])/length(uni$mdn)
# How many percent of rows haven't deviation
length(non_uni$mdn[non_uni$mdn==0])/length(non_uni$mdn)
length(non_uni$mdn[uni$mdn==0])/length(uni$mdn)
# Formal test
wilcox.test(non_uni$mdn, uni$mdn)

## Is removing of short sequences removed diversity too?
pre_non <- non_unif_summary_rows %>% 
  group_by(ortologus_row) %>% 
  summarise(max_codon=(max(table(start_codone))/sum(table(start_codone))),
            codone_name_pre=names(table(start_codone))[which.max(table(start_codone))])
post_non <- non_unif_summary_rows %>%
  group_by(ortologus_row) %>%
  filter(length>=median(length)) %>% 
  summarise(max_codon=(max(table(start_codone))/sum(table(start_codone))),
            codone_name_post=names(table(start_codone))[which.max(table(start_codone))])
diversity_changing <- left_join(pre_non, post_non, by="ortologus_row") %>% 
  mutate(difference=max_codon.y-max_codon.x) %>% 
  mutate(codone_differs=ifelse(codone_name_pre==codone_name_post, "No", "Yes"))


# Alternative start-codons

## Changing in length after new alignment
# Choosing non-uniform genes
non_unif_summary_rows <- summary_rows %>% 
  filter(uniformity == "non-uniform")

# Counting mdn
non_uni <- non_unif_summary_rows %>% 
  group_by(ortologus_row) %>% 
  summarise(max_length_before=max(length),
            med_length_before=median(length),
            mdn_before=normalised_median_deviation(length),
            max_length_after=max(new_length),
            med_length_after=median(new_length),
            mdn_after=normalised_median_deviation(new_length),
            mdn_change=mdn_before-mdn_after)
View(non_uni)

# Comparing maximal frequent start-codons per row before and after length increasing, and after correcting
pre_post_corr_non <- non_unif_summary_rows %>% 
  group_by(ortologus_row) %>% 
  summarise(max_codon_pre=(max(table(start_codone))/sum(table(start_codone))),
            codone_name_pre=names(table(start_codone))[which.max(table(start_codone))],
            max_codon_post=(max(table(new_sc))/sum(table(new_sc))),
            codone_name_post=names(table(new_sc))[which.max(table(new_sc))],
            max_codon_corr=(max(table(new_sc_c))/sum(table(new_sc_c))),
            codone_name_corr=names(table(new_sc_c))[which.max(table(new_sc_c))],
            change_pre_post=max_codon_post-max_codon_pre,
            change_pre_corr=max_codon_corr-max_codon_pre,
            is_diff_pre_post=(codone_name_pre!=codone_name_post),
            is_diff_pre_corr=(codone_name_pre!=codone_name_corr))
View(pre_post_corr_non)

ggplot(summary_rows)+
  geom_histogram(aes(x=delta), binwidth = 50)

# Gene group proportion visualisation
prop_gene_group <- start_codons2 %>% 
  group_by(gene_group, start_type) %>% 
  summarise(count = n()) %>% 
  filter(gene_group != "NA" & is.na(start_type)==FALSE)


positions = c("core", "shell", "cloud")
ggplot(prop_gene_group, aes(x = gene_group, y = count, fill = start_type)) +
  geom_bar(stat = "identity", position="fill")+
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

half_blood_gene <- summary_rows %>%
  filter(ortologus_row==297)
half_blood_gene$product
