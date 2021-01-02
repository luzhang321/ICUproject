# for ICU ITS 
# sig.diff species using wilcox.test;sig.species boxplot in icu+ vs healthy as well as icu- vs healthy comparison;4 candida species boxplot 
# name : cuiqin 
# aim : sig.diff species from species table (metagenomseq normalized & 10% prevelance) 

library(textshape)
library(data.table)
library(tidyverse)



CVSICUm_w <- list()
CVSICUp_w <- list()
ImVSIp_w <- list()

#setwd("/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/")
setwd("/media/lu/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/")


# READ FILE (row:species;col:samples)
ICU_sel_species_t <- fread("result/dada_species_10per.txt") %>%
  as_tibble() %>%
  column_to_rownames(.,var="V1") %>%
  t(.) 



#ICU pheno 

#ICU_pheno <- fread(file = "/Volumes/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt",header = T)
ICU_pheno <- fread(file = "/media/lu//Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt",header = T)



ICU_pheno$group <- factor(ICU_pheno$group, levels = c("Control", "ICU-", "ICU+Mero","ICU+Pip"))
ICU_pheno$Bgroup <- factor(ICU_pheno$Bgroup, levels = c("Control", "ICU-", "ICU+"))

ICU_pheno <- ICU_pheno[-c(74),] # remove AT081



# combine phenotype 
ICU_sel_species_t_groups <- cbind(ICU_sel_species_t[sort(rownames(ICU_sel_species_t)),], ICU_pheno)
ICU_sel_species_t_groups$Bgroup <- as.character(ICU_sel_species_t_groups$Bgroup)
ICU_sel_species_t_groups$Group_num <- ICU_sel_species_t_groups$Bgroup

species_table <- 
  ICU_sel_species_t_groups %>% 
  as_tibble () %>% 
  dplyr::select (Group_num, everything(.))


species_table$Group_num


# kruskal wallis : even though i wrote it here, but i didn't use it, i used separate wilcox.test result for sig.list separately
kruskal_pvalue <- list()
for (i in colnames(species_table)[2:110]){
  p <- kruskal.test(species_table[[i]] ~ Bgroup, data = species_table)$p.value
  kruskal_pvalue[i] <- p
}

kruskal_wallis <- as.data.frame(do.call(rbind,kruskal_pvalue))


CVSICUm_w <- list()
CVSICUp_w <- list()
ImVSIp_w <- list()


for (i in colnames(species_table)[2:(ncol(species_table)-3)]){
  # pairwise.wilcox.test(species_table[[2]], species_table[[1]])
  #        Control ICU-
  #ICU- 1.0000000   NA
  #ICU+ 0.9764976    1
  pvalue_all <- pairwise.wilcox.test(species_table[[i]], p.adjust.method = "none", species_table[[1]])[[3]][1:4]
  #[1] 1.0000000 0.9764976        NA 1.0000000
  #(ICU- vs Control; ICU+ vs Control; ICU- vs ICU+)
  CVSICUm_w[i] <- pvalue_all[1]
  CVSICUp_w[i] <- pvalue_all[2]
  ImVSIp_w[i] <- pvalue_all[4]
}

DF_1.wil <- do.call(rbind, lapply(CVSICUm_w, data.frame))
DF_2.wil <- do.call(rbind, lapply(CVSICUp_w, data.frame))
DF_3.wil <- do.call(rbind, lapply(ImVSIp_w, data.frame))
DF.wil <- cbind(DF_1.wil, DF_2.wil, DF_3.wil)
DF.wil<- na.omit(DF.wil) # remove from 109 to 109
colnames(DF.wil) <- c("CVSM","CVSP","MVSP")
DF.wil$species <- rownames(DF.wil)

# select significant species 
#==========================================

DF.wil$CVSM_fdr <- p.adjust(DF.wil$CVSM, method = "fdr")
DF.wil$CVSP_fdr <- p.adjust(DF.wil$CVSP, method = "fdr")
DF.wil$MVSP_fdr <- p.adjust(DF.wil$MVSP, method = "fdr")


DF.table.WIL <- as.data.table(DF.wil)
# For ICU+ VS ICU- 
siginicu.WIL <- DF.table.WIL[MVSP <= 0.05]$species # 3  length(siginicu.WIL)
# For ICU- VS Control  
sigincVSm.WIL <- DF.table.WIL[CVSM <= 0.05]$species # 11 length(sigincVSm.WIL)
# For ICU+ VS Control
sigincVSp.WIL <- DF.table.WIL[CVSP <= 0.05]$species # 19 length(sigincVSp.WIL)

DF.table.WIL$S1 <- "No"
DF.table.WIL$S2 <- "No"
DF.table.WIL$S3 <- "No"
DF.table.WIL[DF.table.WIL$species %in% siginicu.WIL, ]$S1 <- "Yes"    
DF.table.WIL[DF.table.WIL$species %in% sigincVSp.WIL, ]$S2 <- "Yes"    
DF.table.WIL[DF.table.WIL$species %in% sigincVSm.WIL, ]$S3 <- "Yes"  
nrow(DF.table.WIL[S1=="No" & S2=="No" & S3 == "No"]) #85
DF.table.WIL <- DF.table.WIL[!(S1=="No" & S2=="No" & S3 == "No")] #24 

write_delim(DF.table.WIL, "result/sig_species.txt", delim = "\t")
#DF.table.WIL$species <- gsub("NA", "_sp", DF.table.WIL$species)
#kruskal_wallis$name <- rownames(kruskal_wallis)

#kruskal_wallis$fdr <- p.adjust(kruskal_wallis$V1, method = "fdr")
#ks_value <- kruskal_wallis[DF.table.WIL$species,]



cols <- ncol(species_table)-2

species_df <- species_table %>% 
  dplyr::select (sample, everything ()) %>% 
  as.data.frame(.) %>%
  .[,c(1,3:cols)] # remove Group num column 

rownames(species_df) <- species_df$sample
species_df <- t(species_df[,2:ncol(species_df)]) # remove sample column and transpose 

species_df.WIL <- species_df[rownames(species_df) %in% DF.table.WIL$species,]#(only extract the species that are significantly different): abundance table 

# ### after change the name it will not situable for later analysis 
table(c(rownames(species_df.WIL),DF.table.WIL$species))

rownames(species_df.WIL) <- gsub("NA"," sp",rownames(species_df.WIL))
rownames(species_df.WIL) <- gsub("_"," ",rownames(species_df.WIL))

species_df.WIL.df <- as.data.frame(species_df.WIL) %>%
  rownames_to_column(., var = "species")
  
write_delim(species_df.WIL.df, "result/significant_wilcoxon_table.txt")


# =============select the species that are significantly increased in ICU+ group 
# for example : wilcoxon test 

ICU_sel_sig_species <- list()

ICU_sel_species_table <- as.data.table(t(ICU_sel_species_t))
ICU_sel_species_table$species <- rownames(t(ICU_sel_species_t))


library(ggpubr)


# compare by using median -------------------------------------------------------------
#--------------------------------------------------------------------------------------
# compare the abundance : through the median 
library(ggpubr)


p_ICUp_Control <- list()
p_ICUp_ICUm <- list()
DF.table.WIL.export <- as.data.table(DF.wil) # all species
DF.table.WIL.export$Control <- -0.1
DF.table.WIL.export$ICUp <- -0.1
DF.table.WIL.export$ICUm <- -0.1
DF.table.WIL.export$ICUp_Control <- "NA" # Compare ICUplus and Control, mark yes if higher ICU+
DF.table.WIL.export$ICUp_ICUm <- "NA" # compare ICUplus and ICUminus, mark yes if higher ICU+
DF.table.WIL.export$ICUm_Control <- "NA"  # compare ICUminus and Control 
DF.table.WIL.export$ICUp_lower_Control <- "NA" # compare ICUplus and Control,mark no, if higher in Control
DF.table.WIL.export$ICUp_lower_ICUm <- "NA" # compare ICUplus and ICUminus,mark no, if higher in ICU-


# species_df.WIL : is only the significant species 
# ICU_sel_species_table : abundance table for species 

# compare the abundance : through the median 
for (i in rownames(species_df.WIL)){
  print(i)
  ICU_sel_sig_species[[i]] <- ICU_sel_species_table[species == i]
  test <- ICU_sel_sig_species[[i]]
  #test <- ICU_sel_species_table[species == i] 
  test <- data.frame(test)
  rownames(test) <- test$species
  test <- data.frame(t(test[,1:(ncol(test)-1)])) # transponse data
  test$sample <- rownames(test)
  test_df <- merge(test, ICU_pheno, id = sample) # add group 
  # compare the median of the test 
  test_mean <- test_df[,c(2,4)] %>% group_by(Bgroup) %>% summarise(mean_1 = median(.data[[i]])) %>% ungroup() 
  # add the median to the DF.table.WIL.export: even though here is named by mean, actually it is median  
  DF.table.WIL.export[species == i]$Control <- test_mean$mean_1[1]
  DF.table.WIL.export[species == i]$ICUm <- test_mean$mean_1[2]
  DF.table.WIL.export[species == i]$ICUp <- test_mean$mean_1[3]
  
  ICUp_Control <- test_mean$mean_1[3] - test_mean$mean_1[1]
  ICUm_Control <- test_mean$mean_1[2] - test_mean$mean_1[1]
  if (ICUp_Control > 0){
    #my_comparisons <- list(c("Control","ICU-"),c("Control","ICU+"),c("ICU-","ICU+"))
    #colorinfo = c("#BC3C29FF","#0072B5FF","#E18727FF")
    #stat_compare_means(comparisons=my_comparisons,size = 5)+ 
    #stat_compare_means(label.y.npc = "top", label.x.npc = "right", size =5) +  # Add global p-value + 
    DF.table.WIL.export[species == i]$ICUp_Control <- "Yes"
  }
  if (ICUm_Control >0){
    DF.table.WIL.export[species == i]$ICUm_Control <- "Yes"
  }
  if (ICUm_Control < 0){
    DF.table.WIL.export[species == i]$ICUm_Control <- "No"
  }
  if (ICUp_Control < 0){
    DF.table.WIL.export[species == i]$ICUp_Control <- "No"
  }
  ICUmVSICUp <- test_mean$mean_1[3] - test_mean$mean_1[2]
  if (ICUmVSICUp > 0){
    DF.table.WIL.export[species == i]$ICUp_ICUm <- "Yes"
  }
  if (ICUmVSICUp < 0){
    DF.table.WIL.export[species == i]$ICUp_lower_ICUm <- "No"
  }
} 

DF.table.WIL.export <- DF.table.WIL.export[Control>=0] # only extract those with median : 24

ICUp_Control_higher <- DF.table.WIL.export[ICUp_Control == "Yes"]$species %>% 
  intersect(.,sigincVSp.WIL)  # 4 species 
ICUp_Control_lower <- DF.table.WIL.export[ICUp_Control == "No"]$species %>% 
  intersect(.,sigincVSp.WIL) # 10 species 


ICUp_ICUm_higher <- DF.table.WIL.export[ICUp_ICUm == "Yes"]$species %>% 
  intersect(.,siginicu.WIL) # 0 species

ICUp_ICUm_lower <- DF.table.WIL.export[ICUp_lower_ICUm == "No"]$species %>% 
  intersect(.,siginicu.WIL) # 0 species


ICUm_Control_higher <- DF.table.WIL.export[ICUm_Control == "Yes"]$species %>% 
  intersect(.,sigincVSm.WIL) # 2 species 
ICUm_Control_lower <- DF.table.WIL.export[ICUm_Control == "No"]$species %>%
  intersect(.,sigincVSm.WIL) # 6 species 


# build yes and no label

ICUp_Control_h_l <- data.frame(species = union(ICUp_Control_higher, ICUp_Control_lower), Level = c(rep("Yes",length(ICUp_Control_higher)),rep("No",length(ICUp_Control_lower))), stringsAsFactors = F) # 14species 

ICUp_ICUm_h_l <- data.frame(species = union(ICUp_ICUm_higher, ICUp_ICUm_lower), Level = c(rep("Yes", length(ICUp_ICUm_higher)), rep("No",length(ICUp_ICUm_lower))), stringsAsFactors = F) # 0 species 

ICUm_Control_h_l <- data.frame(species = union(ICUm_Control_higher,ICUm_Control_lower),Level = c(rep("Yes",length(ICUm_Control_higher)),rep("No",length(ICUm_Control_lower))), stringsAsFactors = F) # 8 species 

a_mark_species <- ifelse(ICUp_Control_h_l$Level == "Yes","#D55E00","#0072B2") #48cae4
b <- ifelse(ICUp_ICUm_h_l$Level == "Yes", "#D55E00", "#0072B2")
c_mark_species <- ifelse(ICUm_Control_h_l$Level == "Yes","#D55E00","#0072B2")

# select the specific rows from original species table 
box_plot_ICU <- ICU_sel_species_table %>% as.data.frame()
rownames(box_plot_ICU) <- box_plot_ICU$species

box_plot_ICU <- t(box_plot_ICU[,1:(ncol(box_plot_ICU)-1)])
box_plot_ICU <- as.data.frame(box_plot_ICU)  
box_plot_ICU$sample <- rownames(box_plot_ICU)
box_plot_ICU <- merge(box_plot_ICU, ICU_pheno, id = samples, sort = F) %>% data.table()

icu_control <- box_plot_ICU[,ICUp_Control_h_l$species,with= FALSE] %>% setcolorder(.,ICUp_Control_h_l$species) %>% cbind(.,box_plot_ICU$Bgroup)
#https://stackoverflow.com/questions/12232041/how-to-reorder-data-table-columns-without-copying
#https://stackoverflow.com/questions/32184252/how-to-select-columns-in-data-table-using-a-character-vector-of-certain-column-n

icu_control_melt <- melt(icu_control) 
colnames(icu_control_melt) <- c("Groups","Species","Abundance")
icu_control_melt <- icu_control_melt[Groups != "ICU-"] # remove ICU-

icu_icum <- box_plot_ICU[,ICUp_ICUm_h_l$species,with= FALSE] %>% 
  setcolorder(.,ICUp_ICUm_h_l$species) %>%
  cbind(.,box_plot_ICU$Bgroup)

icu_icum_melt <- melt(icu_icum) 
colnames(icu_icum_melt) <- c("Groups","Species","Abundance")
icu_icum_melt <- icu_icum_melt[Groups != "Control"]

icum_control <- box_plot_ICU[,ICUm_Control_h_l$species,with= FALSE] %>% 
  setcolorder(.,ICUm_Control_h_l$species) %>%
  cbind(.,box_plot_ICU$Bgroup)

icum_control_melt <- melt(icum_control) 
colnames(icum_control_melt) <- c("Groups","Species","Abundance")
icum_control_melt <- icum_control_melt[Groups != "ICU+"]

# calculate log10abundance 
icu_control_melt$logabundance <- log10(icu_control_melt$Abundance+min(icu_control_melt$Abundance[icu_control_melt$Abundance > 0]))
icum_control_melt$logabundance <- log10(icum_control_melt$Abundance+min(icum_control_melt$Abundance[icum_control_melt$Abundance > 0]))

# check pvalue 
wilcox.test(icu_control_melt[Species == "Candida"]$Abundance~icu_control_melt[Species == "Candida"]$Groups)

wilcox.test(icu_control[V2 %in% c("Control","ICU+")][[2]]~icu_control[V2 %in% c("Control","ICU+")][[ncol(icu_control)]])



# produce pvalue 
pvalue <- vector()
for (i in 1:(ncol(icu_control)-1)){
  pvalue[i] <- wilcox.test(icu_control[V2 %in% c("Control","ICU+")][[i]]~icu_control[V2 %in% c("Control","ICU+")][[ncol(icu_control)]])$p.value
}


pvalue_icum <- vector()
for (i in 1:(ncol(icum_control)-1)){
  pvalue_icum[i] <- wilcox.test(icum_control[V2 %in% c("Control","ICU-")][[i]]~icum_control[V2 %in% c("Control","ICU-")][[ncol(icum_control)]])$p.value
}



pvalue_mark <- cut(pvalue,breaks=c(-Inf, 0.001, 0.01, 0.05, 1),label=c("*** ", "** ", "*  ","ns"))
pvalue_mark_icum <- cut(pvalue_icum,breaks=c(-Inf, 0.001, 0.01, 0.05, 1),label=c("*** ", "** ", "*  ","ns"))


main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=9),
                   axis.title=element_text(size=10,face="bold"),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=8),
                   text=element_text(family="sans", size=8),
                   plot.title = element_text(hjust = 0.5,face="bold"))

# supplementary fig.5b
icu_control_melt$Species <- gsub("_"," ",icu_control_melt$Species)
icu_control_melt$Species <- gsub("NA"," sp",icu_control_melt$Species)

icu_control_melt[icu_control_melt$Groups == "Control"]$Groups <- "Healthy"

icum_control_melt$Species <- gsub("_"," ",icum_control_melt$Species)
icum_control_melt$Species <- gsub("NA"," sp",icum_control_melt$Species)
icum_control_melt[icum_control_melt$Groups == "Control"]$Groups <- "Healthy"



a_mark_species <- ifelse(ICUp_Control_h_l$Level == "Yes","#f72585","#48cae4")
b <- ifelse(ICUp_ICUm_h_l$Level == "Yes", "#D55E00", "#0072B2")
c_mark_species <- ifelse(ICUm_Control_h_l$Level == "Yes","#f72585","#48cae4")
# c_mark_species <- ifelse(ICUm_Control_h_l$Level == "Yes","#D55E00","#0072B2")

p <- ggboxplot(icu_control_melt, "Species", "logabundance", orientation = "horizontal",
               fill = "Groups", palette = c("#3b57a6","#f6921d")) + 
  theme(legend.position = "bottom") + ylab("Log10(Abundance)") +
  #stat_compare_means(aes(group = Groups), label = "p.signif") +
  #scale_y_continuous(limits = c(0,250), breaks = c(0,1,3,5,50,150,200,250)) + 
  ggtitle("Healthy vs ICU+") + 
  main_theme +
  labs(fill = "Group") + 
  theme(axis.text.y = element_text(colour=a_mark_species,size = 14)) +
  xlab("Species") 


p <- p + annotate("text", x = c(1:length(a_mark_species)), y = 2.2, label = pvalue_mark, size = 5)   + 
  theme(legend.position = "bottom")
p_icup_healthy <- p + font("y.text", face = "italic")

p2 <- ggboxplot(icum_control_melt, "Species", "logabundance", orientation = "horizontal",
                fill = "Groups", palette = c("#ed3224","#f6921d")) + 
  theme(legend.position = "bottom") + ylab("Log10(Abundance)") +
  #stat_compare_means(aes(group = Groups), label = "p.signif") +
  #scale_y_continuous(limits = c(0,250), breaks = c(0,1,3,5,50,150,200,250)) + 
  ggtitle("Healthy vs ICU-") + 
  main_theme +
  labs(fill = "Group") + 
  theme(axis.text.y = element_text(colour=c_mark_species,size = 14)) +
  xlab("Species") 
  #scale_fill_manual(values = alpha(c("#FC4E07","#00AFBB"), .6))


p2 <- p2 + annotate("text", x = c(1:length(c_mark_species)), y = 2.2, label = pvalue_mark_icum, size = 5) + theme(legend.position = "bottom") 

p2 <- p2 + font("y.text", face = "italic") 

# supplementary fig.5b

p_sig <- ggarrange(p_icup_healthy,p2)
#ggsave("/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/result//S4.combine_median.png", p_sig,width = 12,height = 8)
ggsave("/media/lu/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/result//S4.combine_median.png", p_sig,width = 12,height = 8)

# for publications 

#pdf("/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/result/S4.combine_median.pdf",width = 12,height = 8)
pdf("/media/lu//Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/result/S4.combine_median.pdf",width = 12,height = 8)
p_sig
dev.off()

# #--------------------------4 sig.diff Candida species boxplot  ------------------------------------------------------------------

candida_species <- species_table[,grepl( "Candida" , names(species_table))] %>%
  cbind(species_table$Group_num,. )
library(tidyr)
data_long <- gather(candida_species, species, abundance, Candida_albicans:CandidaNA, factor_key=TRUE)

data_long
colnames(data_long)[1] <- "Group"
data_long <- filter(data_long, !species %in% "CandidaNA")

# i want to draw group boxplot 

main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=9),
                   axis.title=element_text(size=10,face="bold"),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=8),
                   text=element_text(family="sans", size=8),
                   plot.title = element_text(hjust = 0.5,face="bold"))

data_long$log <- log10(data_long$abundance + min(data_long[data_long$abundance>0,]$abundance))


library(ggsci)
# candida albicans 
data_long <- as.data.table(data_long)
write.csv2(data_long,"result/candida_4species_input.csv",quote = F)
data_long[Group == "Control"]$Group <- "Healthy"
data_long$Group <- factor(data_long$Group, levels = c("Healthy","ICU-","ICU+"))
a <- dplyr::filter(data_long, species %in% "Candida_albicans")

p1 <- ggboxplot(a, x = "Group", y = "abundance",fill = "Group") 

my_comparisons <- list( c("Healthy", "ICU+"), c("ICU-", "Healthy"), c("ICU+", "ICU-") )
p1_ca <- p1 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +
  main_theme + 
  ggtitle("Candida albicans") + 
  scale_fill_manual(values = c("#f6921d","#ed3224","#3b57a6")) 
  #scale_fill_manual(values = alpha(c("#00AFBB","#FC4E07", "#E7B800"), .6))

p1_ca <- p1_ca + font("title",face = "bold.italic",size = 14) + 
  font("xlab", size = 14) + 
  font("ylab", size = 14) +
  font("xy.text", size = 14) 

# candida Candida_glabrata

a <- dplyr::filter(data_long, species %in% "Candida_glabrata")

p1 <- ggboxplot(a, x = "Group", y = "abundance",fill = "Group") 


p1_cg <- p1 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +
  main_theme + 
  ggtitle("Candida glabrata") + 
  main_theme + 
  #scale_fill_manual(values = c("#00AFBB","#FC4E07", "#E7B800")) +
  #scale_fill_manual(values = alpha(c("#00AFBB","#FC4E07", "#E7B800"), .6))
  scale_fill_manual(values = c("#f6921d","#ed3224","#3b57a6")) 
p1_cg <- p1_cg +  font("title",face = "bold.italic",size = 14) + 
  font("xlab", size = 14) + 
  font("ylab", size = 14) +
  font("xy.text", size = 14) 

# candida tropicals 
Candida_tropicalis 

a <- dplyr::filter(data_long, species %in% "Candida_tropicalis")

p1 <- ggboxplot(a, x = "Group", y = "abundance",fill = "Group") 


p1_ct <- p1 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +
  main_theme + 
  ggtitle("Candida tropicalis") + 
  #scale_fill_manual(values = c("#00AFBB","#FC4E07", "#E7B800")) +
  main_theme + 
  scale_fill_manual(values = c("#f6921d","#ed3224","#3b57a6")) 
  #scale_fill_manual(values = alpha(c("#00AFBB","#FC4E07", "#E7B800"), .6))
p1_ct <- p1_ct +  font("title",face = "bold.italic",size = 14) + 
  font("xlab", size = 14) + 
  font("ylab", size = 14) +
  font("xy.text", size = 14)  


# candida pseudolambica
Candida_pseudolambica

a <- dplyr::filter(data_long, species %in% "Candida_pseudolambica")

p1 <- ggboxplot(a, x = "Group", y = "abundance",fill = "Group") 


p1_cp <- p1 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +
  main_theme + 
  ggtitle("Candida pseudolambica") + 
  main_theme + 
  scale_fill_manual(values = c("#f6921d","#ed3224","#3b57a6")) 
  #scale_fill_manual(values = c("#00AFBB","#FC4E07", "#E7B800")) +
  #scale_color_npg() +
  #scale_fill_npg() +
  #scale_fill_manual(values = alpha(c("#00AFBB","#FC4E07", "#E7B800"), .6))

p1_cp <- p1_cp +  font("title",face = "bold.italic",size = 14) + 
  font("xlab", size = 14) + 
  font("ylab", size = 14) +
  font("xy.text", size = 14)  


# ggarrange 

#p_candida_species <- ggarrange(p1_ca, p1_cg,p1_ct,p1_cp,ncol = 1,common.legend = TRUE,legend = FALSE)
#p_candida_species_bar_all <- ggarrange(p_bar_mean, p_candida_species,ncol = 2,common.legend = TRUE)
#ggsave("result/candida_species.png", p_candida_species_bar_all, width = 12, height = 12)

p_candida_species <- ggarrange(p1_ca, p1_cg,p1_ct,p1_cp,common.legend = TRUE,legend = FALSE,ncol = 4)
ggsave("result/candida_species_4.png", p_candida_species, width = 12, height = 3)
pdf("result/candida_species_4.pdf", width = 12, height = 4)
p_candida_species
dev.off()

save.image("result/202007.Significantspecies.RData")
load("result/202007.Significantspecies.RData")
