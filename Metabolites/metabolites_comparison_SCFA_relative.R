# Metabolites
# metabolites comparsion - SCFA  
# aim : using kruskal wallis to get sig.SCFA list & and do wilcox.test in each pair & doing FDR correction
# name : cuiqin 

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(data.table)
library(readxl)

#setwd("/Volumes/Lucy/documents/201910_ICU_ITS/04_Metabolics data/Batch2_metabolites_another_company/")

setwd("/media/lu/Lucy/documents/201910_ICU_ITS/04_Metabolics data/Batch2_metabolites_another_company/")
# SCFA 
#----------------------------------------------------------------------------------------


scfa_quantified <- read_excel(path = "data/SCFA_relative_value.xlsx") # check name and input file with the original scfa_results.xlsx;correct! don't need to check again.
pheno <- read.delim("/media/lu/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt")

#pheno <- read.delim("/Volumes/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt")

# i want to replace all < LOD to 0 (while in relative table, there is no LOD)


replace_lod <- function(x){
  a <- str_replace(x, "<LOD", "0")
  return(a)
}

scfa_quantified_replace <- apply(scfa_quantified, 2, replace_lod) %>%
  as_tibble() %>%
  # combine with the group infor
  inner_join(., pheno) 
  

# replace the name space with "_" 
colnames(scfa_quantified_replace) <- str_replace(colnames(scfa_quantified_replace), " ", "_")


# convert the character to numeric 
scfa_quantified_replace <- apply(scfa_quantified_replace[,2:11], 2, as.numeric) %>%
  as_tibble() %>%
  cbind(.,scfa_quantified_replace[,c(1,12,13)])


# test : if it fits normal distribution => don't fit & can't use t-test 
shapiro_1 <- tapply(scfa_quantified_replace$Acetic_acid, scfa_quantified_replace$Bgroup, shapiro.test)
shapiro_2 <- tapply(scfa_quantified_replace$Formic_acid, scfa_quantified_replace$Bgroup, shapiro.test)
shapiro_3 <- tapply(scfa_quantified_replace$Propionic_acid, scfa_quantified_replace$Bgroup, shapiro.test)
shapiro_4 <- tapply(scfa_quantified_replace$Isobutyric_acid, scfa_quantified_replace$Bgroup, shapiro.test)
shapiro_5 <- tapply(scfa_quantified_replace$Butyric_acid, scfa_quantified_replace$Bgroup, shapiro.test)
shapiro_6 <- tapply(scfa_quantified_replace$Isovaleric_acid, scfa_quantified_replace$Bgroup, shapiro.test)
shapiro_7 <- tapply(scfa_quantified_replace$Valeric_acid, scfa_quantified_replace$Bgroup, shapiro.test)

shapiro_8 <- tapply(scfa_quantified_replace$isocaproic_acid, scfa_quantified_replace$Bgroup, shapiro.test)# all of them are identitical  = 0

shapiro_9 <- tapply(scfa_quantified_replace$Caproic_acid,scfa_quantified_replace$Bgroup, shapiro.test)

shapiro_10 <- tapply(scfa_quantified_replace$Heptanoic_acid,scfa_quantified_replace$Bgroup, shapiro.test) #Error in FUN(X[[i]], ...) : all 'x' values are identical

# statistical test 

kruskal_pvalue <- list()
wilcox_pvalue1 <- list()
wilcox_pvalue2 <- list()
wilcox_pvalue3 <- list()


# kruskal wallis 


for (i in colnames(scfa_quantified_replace)[1:10]){
  p <- kruskal.test(scfa_quantified_replace[,i] ~ Bgroup, data = scfa_quantified_replace)$p.value
  kruskal_pvalue[i] <- p
}

pvalue_summary <- do.call(rbind, kruskal_pvalue)

for (i in colnames(scfa_quantified_replace)[1:10]){
  scfa_quantified_replace[,i] <- as.numeric(scfa_quantified_replace[,i])
  p1 <- wilcox.test(scfa_quantified_replace[,i] ~ Bgroup, data = scfa_quantified_replace, subset = Bgroup %in% c("Control","ICU-"))$p.value
  p2 <- wilcox.test(scfa_quantified_replace[,i] ~ Bgroup, data = scfa_quantified_replace, subset = Bgroup %in% c("Control","ICU+"))$p.value
  p3 <- wilcox.test(scfa_quantified_replace[,i] ~ Bgroup, data = scfa_quantified_replace, subset = Bgroup %in% c("ICU-","ICU+"))$p.value
  wilcox_pvalue1[i] <- p1
  wilcox_pvalue2[i] <- p2
  wilcox_pvalue3[i] <- p3 
}

pvalue_summary <- cbind(pvalue_summary, do.call(rbind, wilcox_pvalue1), do.call(rbind, wilcox_pvalue2), do.call(rbind, wilcox_pvalue3))
colnames(pvalue_summary) <- c("Kruskal_Wallis", "ICU-VSControl", "ICU+VSControl", "ICU-VSICU+")


# fdr adjustment 

ks_value <- pvalue_summary %>% 
  as_tibble() %>% 
  .$Kruskal_Wallis

ks_fdr <- p.adjust(ks_value, method = "fdr")

pvalue_summary_fdr <- pvalue_summary %>%
  as.data.frame() %>%
  mutate(ks_fdr = ks_fdr, scfa = rownames(pvalue_summary))

write_delim(pvalue_summary_fdr,"result/relative/scfa_comparisons.txt",delim = "\t")
write_delim(scfa_quantified_replace,"result/relative/scfa_all.txt",delim = "\t")

# barplot 

SCFA <- reshape2::melt(scfa_quantified_replace[,c(1:10,13)], id = "Bgroup")


colnames(SCFA)[1] <- "Groups"
p_scfa <- ggbarplot(SCFA, x = "variable", y = "value", add = "mean_se",
                    color = "Groups",fill = "Groups", palette = "jco", 
                    position = position_dodge(0.8))+
  stat_compare_means(aes(group = Groups), label = "p.signif") +
  rotate_x_text(angle = 45) + xlab("Short Chain Fatty Acids") + 
  ylab("relative content") +
  theme(legend.position="bottom") #+
  #theme_minimal() 

p_scfa

ggsave("result/relative/scfa_combine_barplot.png", p_scfa, width = 10, height = 6)


