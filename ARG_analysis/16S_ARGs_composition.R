# Resistome Analysis/Profiling 
# ARGs compositon difference among different groups
# name : cuiqin 
# date : 03.25 

# ARG 16S normalized file + pheno file 
# Rscript script.R args group-file outdir outname 

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Rscript ARG_RPKM_difference.R arg-RPKM-file group-file outdir outname ")
} else if (length(args)==4) {
  print(paste("ARG class file:",args[1],";pheon_file,",args[2],",outdir:",args[3],",outname:",args[4],sep = ""))
} 

library(ggplot2)
library(data.table)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(hrbrthemes)
library(stringr)



setwd("/media/lu/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/webserver_identity80_cov50/data/")
#ICU_pheno <- fread("/Volumes/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt", header = T)
#ICU_pheno <- fread("/media/lu/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/tidy_directory/script/Resistome_composition/data/ICU_pheno_anno.txt")
ICU_pheno <- fread(args[2], header = T)
ICU_pheno$Bgroup <- factor(ICU_pheno$Bgroup, levels = c("Control","ICU-","ICU+"))



#ARGs_RPKM <- fread("../data/merged_16S_subtypes_cov50_identity80_website_new.txt",header = T)
ARGs_RPKM <- fread(args[1], header = T)

colnames(ARGs_RPKM) <- str_split(colnames(ARGs_RPKM),pattern = ".deep",simplify = T)[,1]

ARGs_RPKM_df <- dcast(melt(ARGs_RPKM, id =1,variable.name = "sample"), sample ~ ID) %>%
  merge(., ICU_pheno, by = "sample") 


# sum up the values 

ncol <- ncol(ARGs_RPKM_df)
ncol_3 <- ncol - 3

ARGs_RPKM_sum <- ARGs_RPKM_df %>%
  column_to_rownames(., var = "sample") %>%
  mutate(rowsum = rowSums(.[1:ncol_3]))


# test if they fit normal distribution 
tapply(ARGs_RPKM_sum$rowsum, ARGs_RPKM_sum$Bgroup, shapiro.test)

pv_ICUm_ICUp <- wilcox.test(rowsum~Bgroup, data = ARGs_RPKM_sum, subset = Bgroup %in% c("ICU-","ICU+"))[3]$p.value

pv_ICUp_Control <- wilcox.test(rowsum~Bgroup, data = ARGs_RPKM_sum, subset = Bgroup %in% c("Control","ICU+"))[3]$p.value

pv_ICUm_Control <- wilcox.test(rowsum~Bgroup, data = ARGs_RPKM_sum, subset = Bgroup %in% c("Control","ICU-"))[3]$p.value

pv_all <- kruskal.test(rowsum~Bgroup,data= ARGs_RPKM_sum)[3]$p.value


# convert to log value 
ARGs_RPKM_sum$log10value <- log10(ARGs_RPKM_sum$rowsum)
ARGs_RPKM_sum$Group <- ARGs_RPKM_sum$Bgroup 
ARGs_RPKM_sum <- as.data.table(ARGs_RPKM_sum)
ARGs_RPKM_sum[Group == "Control"]$Group <- "Healthy"
ARGs_RPKM_sum$Group <- factor(ARGs_RPKM_sum$Group, levels = c("Healthy","ICU-","ICU+"))



# plot theme 

main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=9),
                   #axis.title=element_text(size=10,face="bold"),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=8),
                   text=element_text(family="sans", size=10),
                   plot.title = element_text(hjust = 0.5,face="bold"))


# boxplot 

plot_box <- function(data_frame,title_name,p_value){
  p <- data_frame %>%
    ggplot( aes(x=Group, y=log10value, fill = Group), show.legend = FALSE) +
    #ggplot( aes(x=Group, y=rowsum, fill = Group), show.legend = FALSE) + 
    geom_boxplot() +
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.2, alpha=0.3,inherit.aes = TRUE) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    #theme_ipsum() +
    #ggtitle(title_name) +
    xlab("")+
    ylab("Abundance (log10(RPKM))") +
    theme(legend.position="none") + 
    main_theme + 
    annotate("text",x = 1.3, y = 3.8, hjust =1, label = paste("p-value:",formatC(p_value,format = "e",digits = 2),sep = ""),size =3) 
    #geom_boxplot(fill='#A4A4A4', color="darkred") 
    return(p) 
}

p_all <- plot_box(ARGs_RPKM_sum, "ICUm VS ICUp VS Control",pv_all) + scale_fill_manual(values = c("#00AFBB","#FC4E07", "#E7B800"),labels = c("Healthy","ICU-","ICU+"),name = "Group")

p_ICUp_ICUm <- ARGs_RPKM_sum %>% .[Group %in% c("ICU-","ICU+")] %>% plot_box(., "ICU- vs ICU+",pv_ICUm_ICUp) + scale_fill_manual(values = c("#FC4E07", "#E7B800"),labels = c("ICU-","ICU+"),name = "Group")

p_ICUp_Control <- ARGs_RPKM_sum %>% .[Group %in% c("ICU+","Healthy")] %>% plot_box(., "Healthy vs ICU+",pv_ICUp_Control) + scale_fill_manual(values = c("#00AFBB", "#E7B800"),labels = c("Healthy","ICU+"),name = "Group") 

p_ICUm_Control <- ARGs_RPKM_sum %>% .[Group %in% c("ICU-","Healthy")] %>% plot_box(., "Healthy vs ICU-",pv_ICUm_Control) + scale_fill_manual(values = c("#00AFBB", "#FC4E07"),labels = c("Healthy","ICU-"),name = "Group")

p_combine_three <- ggarrange(p_ICUp_ICUm,p_ICUp_Control,p_ICUm_Control,ncol = 3)
#ggsave("mean_ARG_length/S2.composition_ARGs.png",p_combine_three,width = 12,height = 4.06)

ggsave(paste(args[3],args[4],".S2.composition_ARGs_boxplot_accumulative.png",sep= ""),p_combine_three, width = 12, height = 4.06)


# violin plot 

plot_violin <- function(data_frame,title_name,p_value){
  p <- data_frame %>%
    ggplot( aes(x=Group, y=log10value, fill = Group), show.legend = FALSE) +
    #ggplot( aes(x=Group, y=rowsum, fill = Group), show.legend = FALSE) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white")+
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    #theme_ipsum() +
    #ggtitle(title_name) +
    xlab("")+
    ylab("log10(16S Normalized Abundance)") +
    theme(legend.position="none") + 
    main_theme + 
    annotate("text",x = 1.3, y = 3.8, hjust =1, label = paste("p-value:",formatC(p_value,format = "e",digits = 2),sep = ""),size =3) 
  #geom_boxplot(fill='#A4A4A4', color="darkred") 
  return(p) 
}


p_ICUp_ICUm_violin <- ARGs_RPKM_sum %>% .[Group %in% c("ICU-","ICU+")] %>% plot_violin(., "ICU- vs ICU+",pv_ICUm_ICUp) + scale_fill_manual(values = c("#ed3224", "#3b57a6"),labels = c("ICU-","ICU+"),name = "Group")

p_ICUp_Control_violin <- ARGs_RPKM_sum %>% .[Group %in% c("ICU+","Healthy")] %>% plot_violin(., "Healthy vs ICU+",pv_ICUp_Control) + scale_fill_manual(values = c("#f6921d", "#3b57a6"),labels = c("Healthy","ICU+"),name = "Group") 

p_ICUm_Control_violin <- ARGs_RPKM_sum %>% .[Group %in% c("ICU-","Healthy")] %>% plot_violin(., "Healthy vs ICU-",pv_ICUm_Control) + scale_fill_manual(values = c("#f6921d", "#ed3224"),labels = c("Healthy","ICU-"),name = "Group")

p_combine_three_violin <- ggarrange(p_ICUp_ICUm_violin,p_ICUm_Control_violin,p_ICUp_Control_violin,ncol = 3)

pdf("../script/202006_modification_plot/S2.composition_ARGs_violinplot.pdf",width = 12, height = 4.06)
p_combine_three_violin
dev.off()

#ggsave(paste("../script/202006_modification_plot/S2.composition_ARGs_violinplot_new.png",sep= ""),p_combine_three_violin, width = 12, height = 4.06)

ggsave(paste(args[3],args[4],"S2.composition_ARGs_violinplot_new.png",sep= ""), p_combine_three_violin, width = 12, height = 4.06)









