# Resistome Analysis/Profiling 
# ARGs compositon difference among different groups: beta-diversity bray-curtis adnois 
# name : cuiqin date : 03.25 (add parameters 04.29)

library(ggplot2)
library(vegan) # version 2.5-6
library(magrittr)
library(data.table)
library(stringr)
library(ggpubr)
library(tidyverse)
library(plyr)

# Rscript script.R args-file group-file outdir outname 

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Rscript beta-diversity.R arg-RPKM-file group-file outdir outname ")
} else if (length(args)==4) {
  print(paste("ARG family file:",args[1],";pheon_file,",args[2],",outdir:",args[3],",outname:",args[4],sep = ""))
} 


# ARG file + pheno file
#setwd("/Volumes/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/identity_80_deeparg/output/")
setwd("/media/lu/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/webserver_identity80_cov50/")
#ICU_pheno <- fread("/Volumes/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt", header = T)
ICU_pheno <- fread("/media/lu/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt", header = T)

ICU_pheno <- fread(args[2], header = T)

ICU_pheno$Bgroup <- factor(ICU_pheno$Bgroup, levels = c("Control","ICU-","ICU+"))



ARGs_RPKM <- fread("data/merged_16S_subtypes_cov50_identity80_website_new.txt", header = T)
ARGs_RPKM <- fread(args[1])



colnames(ARGs_RPKM) <- str_split(colnames(ARGs_RPKM),pattern = ".deep",simplify = T)[,1]

ARGs_RPKM_tmp <- dcast(melt(ARGs_RPKM, id =1,variable.name = "sample"), sample ~ ID) %>% as.data.frame()

write_delim(ARGs_RPKM_tmp,paste(args[3],args[4],"RPKM_transponse.txt",sep = ""),delim = "\t")

rownames(ARGs_RPKM_tmp) <- ARGs_RPKM_tmp[,1]
ARGs_RPKM_tmp <- ARGs_RPKM_tmp[,-1]




# plot 

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



PcoA_plots <- function(points_df, matrix_t){
  adonis_result <- adonis(matrix_t ~ Bgroup, data = points_df, distance = "bray")
  print(adonis_result)
  pvalue <- adonis_result$aov.tab$`Pr(>F)`[1]
  F_value <- round(adonis_result$aov.tab$F.Model[1], 3)
  p <- ggplot(points_df, aes(x=x, y=y, color=Bgroup))
  p <- p + geom_point(alpha=.7, size=2) +
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
         color = "Group") + main_theme + 
    #stat_ellipse(type = "norm") +  # draw a circle 
    #geom_polygon(data = points_df,alpha = 0.2) +
    annotate(geom = "text", label = paste("F:",F_value, "; pvalue:",pvalue, sep = " "),  x=max(points_df$x), y = min(points_df$y)-0.05,hjust=1,  color = "black", size = 3)
  q <- p #scale_color_brewer(palette = "Set1")
  return(q)
}


# bray-curits + pcoa axis calculation 
set.seed(9)
ARGs_RPKM_bray_curtis <- vegdist(ARGs_RPKM_tmp, method = "bray", diag = T, upper = T)
set.seed(4)
pcoa_type <- cmdscale(ARGs_RPKM_bray_curtis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points_type <- as.data.frame(pcoa_type$points) # get coordinate string, format to dataframme
colnames(points_type) = c("x", "y", "z") 
eig <- pcoa_type$eig
# add sample infomation 
points_type <- points_type %>% 
  data.frame(sample = rownames(.)) %>%
  merge(., ICU_pheno, id = sample, sort = FALSE)


Beta_all <- PcoA_plots(points_type, ARGs_RPKM_tmp) + scale_color_manual(values = c("#00AFBB","#FC4E07", "#E7B800"),labels = c("Healthy","ICU-","ICU+"),name = "Group")

# two group comparison 
two_groups_pcoa <- function(points_type,matrix,a,b){
  matrix$sample <- rownames(matrix)
  matrix <- as.data.table(matrix)
  points_sel <- points_type %>% 
    as.data.table(.) %>% 
    .[Bgroup %in% c(a,b)]
  matrix_sel <- matrix[sample %in% points_sel$sample] %>% as.data.frame()
  rownames(matrix_sel) <- matrix_sel$sample
  matrix_sel <- matrix_sel[,-ncol(matrix_sel)]
  
  df <- points_sel
  find_hull <- function(df) df[chull(df$x, df$y),]
  hulls <- ddply(points_sel, "Bgroup", find_hull)
  
  Beta <- PcoA_plots(points_sel, matrix_sel) + geom_polygon(data = hulls,alpha = 0.2)
    
  return(Beta)
}

set.seed(4)
ICUplus_control <- two_groups_pcoa(points_type, ARGs_RPKM_tmp, "ICU+","Control") + scale_color_manual(values = c("#f6921d", "#3b57a6"),labels = c("Healthy","ICU+"),name = "Group")+main_theme 

set.seed(4)
ICUminus_control <- two_groups_pcoa(points_type, ARGs_RPKM_tmp, "ICU-","Control") + scale_color_manual(values = c("#f6921d", "#ed3224"),labels = c("Healthy","ICU-"),name = "Group") + main_theme
set.seed(4)
ICUminus_ICUplus <- two_groups_pcoa(points_type, ARGs_RPKM_tmp, "ICU+","ICU-") + scale_color_manual(values = c("#ed3224", "#3b57a6"),labels = c("ICU-","ICU+"),name = "Group") + main_theme

Beta_diversity <- ggarrange(Beta_all, ICUplus_control, ICUminus_control,ICUminus_ICUplus)
#ggsave(Beta_diversity, filename = "mean_ARG_length/S4.Betadiveristy.Bray-curtis.png")
Beta_diversity <- ggarrange(ICUplus_control, ICUminus_control,ICUminus_ICUplus,ncol = 3)

ggsave(Beta_diversity, filename = paste(args[3],args[4],".S3.beta_diversity.Bray-curtis.png",sep=""))

pdf("script/202006_modification_plot/S3.beta_diversity.bray-curtis.pdf",width = 12, height = 4.06)# newest 
Beta_diversity
dev.off()
