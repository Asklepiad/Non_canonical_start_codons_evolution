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

options(scipen = 999)

# Boxplot analyzer
data_joiner <- function(path){
  setwd(path)
  df <- list.files(path=path, pattern = "*_for_common_boxplot.csv", full.names = TRUE) %>% 
    lapply(read_csv) %>%
    do.call(rbind, .)
  dir.create(path=file.path(path, "combined"))
  write_csv(df, file=file.path(path, "combined","combined.csv"))
}

path = "/home/asklepiad/bioinf/start_codons/BI_project_2022/C_psittaci_folder/boxplots/"
data_joiner(path)
boxplots_comm <- read_csv(file.path(path, "combined","combined.csv"))
boxplots_comm$start_codon <- as.factor(boxplots_comm$start_codon)
boxplots_comm$specie <- as.factor(boxplots_comm$specie)

order <- c("S_ruber", "M_bovis", "C_pseudotub", "B_bifidum", "B_thuringie", 
           "B_subtilis", "B_pumilus", "S_simulans", "S_haemolyti", "L_lactis", 
           "E_faecalis", "W_cibaria", "L_salivariu", "C_trachomat", "B_cenocepac", 
           "N_meningiti", "N_gonorrhoe", "A_tumefacie", "H_somni", "V_anguillar",  
           "V_campbelli", "V_cholerae", "E_hormaeche", "K_pneumonia", "E_coli", 
           "A_baumannii", "S_maltophil", "X_fastidios",  "X_campestri", 
          "P_putida", "P_fluoresce", "P_aeruginos", "F_tularensi", "L_interroga")
ggplot(boxplots_comm)+    # Phylogenetic boxplot
  geom_boxplot(aes(x=specie, y=number, color=start_codon))+
  scale_x_discrete(limits=order)+
  theme(axis.text.x = element_text(angle = 35))

sum_of_ass <- boxplots_comm %>% 
  group_by(specie, p_c_unity) %>% 
  summarise(ass_sum=sum(number)) %>% 
  ungroup(specie, p_c_unity)

barplots_comm <- left_join(boxplots_comm, sum_of_ass) %>% 
  mutate(number_rel = number/ass_sum) %>% 
  group_by(specie, start_codon) %>% 
  summarise(mean_percent = mean(number_rel))

ggplot(subset(barplots_comm, start_codon!="ATG"))+     #   Phylogenetic barplot
  geom_col(aes(x=specie, y=mean_percent, fill=start_codon), color="black")+
  scale_x_discrete(limits=order)+
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"),
        legend.position = "none",
        plot.margin = margin(,,2,, "cm"),
        axis.title.x = element_text(vjust=-5))+
  scale_fill_manual(values = c("GTG" = "dark green",
                               "TTG" = "blue"))+
  xlab(label = "Specie")+
  ylab(label = "Proportion of start-codon")+
  labs("fill" = "Legend")
ggsave("../data/barplot_common_2scs.png",  width = 40, height = 20, units = "cm", dpi = 600)

mean_of_org <- sum_of_ass %>% 
  group_by(specie) %>% 
  summarise(org_mean=mean(ass_sum)) %>% 
  arrange(., desc(org_mean))



ggplot(subset(barplots_comm, start_codon!="ATG"))+     #   Genome size barplot
  geom_col(aes(x=specie, y=mean_percent, fill=start_codon), color="black")+
  scale_x_discrete(limits=mean_of_org$specie)+
  theme(axis.text.x = element_text(angle = 35),
        axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))+
  scale_fill_manual(values = c("GTG" = "dark green",
                                "TTG" = "blue"))+
  xlab(label = "Specie")+
  ylab(label = "Proportion of start-codon")+
  labs("fill" = "Legend")
ggsave("../data/barplot_common_2scs.png",  width = 40, height = 20, units = "cm", dpi = 600)

# Cog common stats ####

start_codons2_joiner <- function(path){
  setwd(path)
  df <- list.files(path=path, pattern = "*.csv", full.names = TRUE) %>% 
    lapply(read_csv) %>%
    do.call(rbind, .)
  dir.create(path=file.path(path, "combined"))
  write_csv(df, file=file.path(path, "combined","combined.csv"))
}

path="/home/asklepiad/bioinf/start_codons/BI_project_2022/C_psittaci_folder/orto_all"
start_codons2_joiner(path)
start_codons2_comm <- read_csv(file.path(path, "combined","combined.csv"))
start_codons2_comm$scs_perrow <- start_codons2_comm %>% 
  select(ATG:TTG) %>% 
  apply(., 1, function(x) ifelse(x[1] >= x[2] & x[1] >= x[3], "ATG",
                                 ifelse(x[2] >= x[3] & x[2] > x[1], "GTG", "TTG")))

start_codons2_comm %>% 
  select(cog, scs_perrow, organism) %>% 
  group_by(cog, organism) %>% 
  summarise(ATG=sum(scs_perrow=="ATG"),
            GTG=sum(scs_perrow=="GTG"),
            TTG=sum(scs_perrow=="TTG")) %>% 
  filter(cog %in% mostly_noncanonical$cog)
  
mostly_noncanonical <- start_codons2_comm %>% 
  select(cog, scs_perrow, organism) %>% 
  group_by(cog) %>% 
  summarise(ATG=sum(scs_perrow=="ATG"),
            GTG=sum(scs_perrow=="GTG"),
            TTG=sum(scs_perrow=="TTG")) %>% 
  mutate(atg_ratio=ATG/(ATG+GTG+TTG),
         total=ATG+GTG+TTG) %>% 
  filter(atg_ratio<0.5)

