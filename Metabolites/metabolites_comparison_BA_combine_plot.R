# Metabolites
# Bile Acids comparison + produce combine plot (bile acids & scfa)
# name : cuiqin 
# this one sum up all unknown bile acids & do kruskal wallis & do wilcox.test 


library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(data.table)
library(readxl)
library(readr)

#setwd("/Volumes/Lucy/documents/201910_ICU_ITS/04_Metabolics data/Batch2_metabolites_another_company/")

setwd("/media/lu/Lucy/documents/201910_ICU_ITS/04_Metabolics data/Batch2_metabolites_another_company/")


# Bile Acids 
#----------------------------------------------------------------------------------------

ba_quantified <- read_excel(path = "data/BA_quantified_value.xlsx") # i have checked, it is totally right. 
colnames(ba_quantified)[1] <- "sample"

# i want to join it with the pheno file 

pheno <- read_delim("/media/lu//Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt",delim = "\t")

#pheno <- read_delim("/Volumes/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt",delim = "\t")


ba_quantified <- ba_quantified %>%
  as_tibble() %>%
  # combine with the group infor
  inner_join(., pheno) 


# replace the name space with "_" 
colnames(ba_quantified) <- str_replace(colnames(ba_quantified), " ", "_")


# convert the character to numeric 
ba_quantified <- apply(ba_quantified[,2:52], 2, as.numeric) %>%
  as_tibble() %>%
  cbind(.,ba_quantified[,c(1,53,54)])


# test : if it fits normal distribution 

sharpiro_test <- list()
for (i in 1:51){
  sharpiro_test[[i]] <- tapply(ba_quantified[[i]], ba_quantified$Bgroup, shapiro.test)
} #Error in FUN(X[[i]], ...) : all 'x' values are identical, it doesn't matter,check sharpiro_test[[2]],not all >0.1 didn't fit normality


# sum the like bile acids together 


unknown_ba_sum <- dplyr::select(ba_quantified, contains(c("like","sample","Bgroup"))) %>% # only select like bile acids : 25 
  reshape2::melt(., id.vars = c("sample","Bgroup")) %>% # long format, easily for sum up 
  dplyr::group_by(sample) %>%  # group by sample, prepare for next step : sum up 
  dplyr::summarise(sum = sum(value)) %>%  # sum up by group 
  tibble::add_column(metabolites = 'unknown bile acids') %>% # give a new name unknown bile acids
  tidyr::spread(.,metabolites,sum)

# produce a new bile acid ba_quantified_unknown_sum table 

ba_quantified_unknown_sum <- cbind(unknown_ba_sum[,2],dplyr::select(ba_quantified, !contains("like")))



# statistical test 

kruskal_pvalue <- list()
wilcox_pvalue1 <- list()
wilcox_pvalue2 <- list()
wilcox_pvalue3 <- list()


for (i in colnames(ba_quantified_unknown_sum)[1:27]){ # check colnames(ba_quantified_unknown_sum); the 1-27 is the metabolites column
  p <- kruskal.test(ba_quantified_unknown_sum[,i] ~ Bgroup, data = ba_quantified)$p.value
  kruskal_pvalue[i] <- p
}

pvalue_summary <- do.call(rbind, kruskal_pvalue)

for (i in colnames(ba_quantified_unknown_sum)[1:27]){
  ba_quantified_unknown_sum[,i] <- as.numeric(ba_quantified_unknown_sum[,i])
  p1 <- wilcox.test(ba_quantified_unknown_sum[,i] ~ Bgroup, data = ba_quantified_unknown_sum, subset = Bgroup %in% c("Control","ICU-"))$p.value
  p2 <- wilcox.test(ba_quantified_unknown_sum[,i] ~ Bgroup, data = ba_quantified_unknown_sum, subset = Bgroup %in% c("Control","ICU+"))$p.value
  p3 <- wilcox.test(ba_quantified_unknown_sum[,i] ~ Bgroup, data = ba_quantified_unknown_sum, subset = Bgroup %in% c("ICU-","ICU+"))$p.value
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

pvalue_summary_df <- pvalue_summary %>%
  as.data.frame(.) %>%
  rownames_to_column(., var = "bile_acid") %>%
  mutate(fdr = ks_fdr)

write_delim(pvalue_summary_df,"result/relative_sum_up_all_unknown_bileacids/ba_comparisons.txt",delim = "\t")
pvalue_summary_df_sel <- filter(pvalue_summary_df, fdr < 0.05)
write_delim(pvalue_summary_df_sel,"result/relative_sum_up_all_unknown_bileacids/ba_comparisons_sig.txt",delim = "\t")
write_delim(ba_quantified_unknown_sum,"result/relative_sum_up_all_unknown_bileacids/ba_all.txt",delim = "\t")


# barplot 
# ----------------------------------------------------------------------------
# long format 
ba_quantified_long <- reshape::melt(ba_quantified_unknown_sum[,c(1:27,30)], id.vars = "Bgroup") # cause after sum up, there is only 27 bile acids 

ba_quantified_long_sel <- ba_quantified_long %>%
  dplyr::filter(.,variable %in% pvalue_summary_df_sel$bile_acid)

colnames(ba_quantified_long_sel)[1] <- "Groups"


#some value lower than 0 -ask them already, suggested directly use 

p_ba <- ggbarplot(ba_quantified_long_sel, x = "variable", y = "value", add = "mean_se",
                  color = "Groups",fill = "Groups", palette = "jco", 
                  position = position_dodge(0.8))+
  stat_compare_means(aes(group = Groups), label = "p.signif") +
  rotate_x_text(angle = 45) + xlab("Bile acids") + 
  ylab("content(ug/g)") +
  theme(legend.position="bottom") #+
#theme_minimal() 

ba_quantified_long_sel$log <- log10(ba_quantified_long_sel$value + 1)

p_ba_log <- ggbarplot(ba_quantified_long_sel, x = "variable", y = "log", add = "mean_se",
                    color = "Groups",fill = "Groups", palette = "jco", 
                    position = position_dodge(0.8))+
  stat_compare_means(aes(group = Groups), label = "p.signif") +
  rotate_x_text(angle = 45) + xlab("Bile acids") + 
  ylab("log10(content(ug/g)+1)") +
  theme(legend.position="bottom") #+
#theme_minimal() 

p_ba_log
p_ba

ggsave("result/relative_sum_up_all_unknown_bileacids/p_ba_barplot.png", p_ba, width = 10, height = 26)


ggsave("result/relative_sum_up_all_unknown_bileacids/p_ba_barplot_log.png", p_ba_log, width = 10, height = 6)





# combine them together for figure 5 
# -------------------------------------------------------------------------------------------------------------------------------

scfa <- fread("result/relative/scfa_all.txt")
SCFA <- reshape2::melt(scfa[,c(1:10,13)], id = "Bgroup")
colnames(SCFA) <- c("Groups","variable","value")

SCFA_quantified_long_sel <- filter(SCFA, !variable %in% c("Formic_acid","isocaproic_acid","Isovaleric_acid","Isobutyric_acid")) # Remove un-Sig.Different scfa & un consistent scfa#

SCFA_quantified_long_sel$variable <- str_replace_all(SCFA_quantified_long_sel$variable,"_"," ")

ba_quantified_long_sel$log <- log(ba_quantified_long_sel$value+1)

ba_quantified_long_sel_log <- ba_quantified_long_sel[,c("Groups","variable","log")]
colnames(ba_quantified_long_sel_log)[3] <- "value"
  
combine_metabolites <- rbind(SCFA_quantified_long_sel,ba_quantified_long_sel_log)
colnames(combine_metabolites)[1] <- "Group"
combine_metabolites$variable <- str_replace_all(combine_metabolites$variable,"_"," ")

combine_metabolites <- as.data.table(combine_metabolites)
combine_metabolites[Group == "Control",]$Group <- "Healthy"
combine_metabolites$Group <- factor(combine_metabolites$Group, levels = c("Healthy","ICU-","ICU+"))

write_csv2(combine_metabolites,path = "result/relative_sum_up_all_unknown_bileacids/combine_barplot_input.txt")


# plot read table again 
#-------------------------------------------------------------------------------
library(data.table)
library(readr)
library(magrittr)
combine_metabolites <- read_csv2("result/relative_sum_up_all_unknown_bileacids/combine_barplot_input.txt")

#1. rename the colname 
colnames(combine_metabolites)[2]<-"metabolites"


#2. statistic 


library(rstatix)
stat.test.metabolties <- combine_metabolites %>%
  group_by(metabolites) %>%
  wilcox_test(value ~ Group,p.adjust.method = "none") #Compute easily statistical tests (t_test() or wilcox_test()) using the rstatix package


# the result using real value and log value are the same, because wilcox using rank not the value itself 
# test <- ba_quantified_long_sel %>%
#   select(Groups, variable, value)
# test_sta <- test %>%
#   group_by(variable) %>%
#   wilcox_test(value ~ Groups,p.adjust.method = "none")


# order it 
stat.test.metabolties$metabolites <- factor(stat.test.metabolties$metabolites, levels = unique(combine_metabolites$metabolites))


stat.test.metabolties <- stat.test.metabolties[order(stat.test.metabolties$metabolites),] %>%
  rstatix::add_xy_position(x = "metabolites", fun = "mean_sd", dodge = 0.8) %>%
  dplyr::mutate(y.position = rep(c(6.5,6.8,7.1),13)) # 6 SCFA AND 7 BA


#3. plot part 
#3.1 basic 
p_combine <- ggbarplot(combine_metabolites, x = "metabolites", y = "value", add = "mean_se",
                color = "Group",fill = "Group", palette = c("#f6921d", "#ed3224", "#3b57a6"), 
                position = position_dodge(0.8))+
  #stat_compare_means(aes(group = Group), label = "p.signif") +
  rotate_x_text(angle = 45) + ylab("Relative content") +
  theme(legend.position="right") + 
  theme(legend.key.size = unit(0.3, "cm"))

# 3.2 add statistic : even thought this is p.adj.signif, i used none for adjustment 

p_combine1 <- p_combine + stat_pvalue_manual(
  stat.test.metabolties, label = "p.adj.signif",tip.length = 0.00,size = 1.88,
  bracket.nudge.y = 0
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))




# add second y-axis 

#p_combine2 <- p_combine1 + scale_y_continuous(  "Relative content",  sec.axis = sec_axis(~ . * 1, name = "Log(Absolute content() + 1)") ) + 
  annotate("segment", x = 1, xend = 6, y = 8, yend = 8, colour = "#F7B22A", size=1, alpha=0.6) + annotate("segment", x = 6.5, xend = 13, y = 8, yend = 8, colour = "#FA1846", size=1, alpha=0.6) 


p_combine2 <- p_combine1 + scale_y_continuous(  "Relative content",  sec.axis = sec_axis(~ . * 1, name = expression("Log(Absolute content" ~ (mu ~ g / g ) ~ "+1 )"))) + 
  annotate("segment", x = 1, xend = 6, y = 8, yend = 8, colour = "#F7B22A", size=1, alpha=0.6) + annotate("segment", x = 6.5, xend = 13, y = 8, yend = 8, colour = "#FA1846", size=1, alpha=0.6) 




annotation <- data.frame(
  x = c(4.5,9.5),
  y = c(8.5,8.5),
  label = c("SCFA", "BA")
)

main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7),
                   plot.title = element_text(hjust = 0.5,face="bold"))
# Add text
p_combine3 <- p_combine2 + geom_text(data=annotation, aes( x=x, y=y, label=label),                 
                     color="black", 
                     size=2 , fontface="bold" ) + xlab("") + main_theme 
  

pdf("result/relative_sum_up_all_unknown_bileacids/Figure5.metabolites_combine.pdf",width = 10, height = 4)
p_combine3
dev.off()



# check which one is lowest 
# ----------------------------------------------------------------------------------------------------
# because in the plot it's easy to see which one ICU+  is higher or lower, except the ones below 

# SCFA except the hepatonotic 

filter(SCFA_quantified_long_sel,variable == "Heptanoic_acid") %>%
  group_by(Groups) %>%
  summarise_at("value",mean)


#   Groups   value
#<fct>    <dbl>
#  1 Control 0.0736
#2 ICU-    0.0217
#3 ICU+    0.0211

#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/79-plot-meansmedians-and-error-bars/
# the figure is shown the mean+/- se , because use the add(mean_se) parameter 
