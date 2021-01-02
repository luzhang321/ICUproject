# for ICU ITS 
# 24 sig.species heatmap 
# name : cuiqin 

# draw heatmap 
library(pheatmap)
library(tidyverse)
library("RColorBrewer")
library(data.table)
library(tibble)

#setwd("/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS")
setwd("/media/lu/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS")

# read input data 

abundance <- fread("result/significant_wilcoxon_table.txt") %>%
  column_to_rownames(., var = "species")


pheatmap(abundance)

# group information 

group <- fread("result/significant species table_benchling.csv")



# pheno type 

#ICU_pheno <- fread(file = "/Volumes/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt",header = T)
ICU_pheno <- fread(file = "/media/lu/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt",header = T)

ICU_pheno <- ICU_pheno[-c(74),] # remove AT081

ICU_pheno$Group <- ICU_pheno$Bgroup
ICU_pheno[Group == "Control"]$Group <- "Healthy"
ICU_pheno$Group <- factor(ICU_pheno$Group, levels = c("Healthy","ICU-","ICU+"))



#######################################################

# complex heatamp
library(ComplexHeatmap)
#devtools::install_github("caleblareau/BuenColors")

library("BuenColors")

col <- jdb_color_maps[1:6]

library(circlize)

library("RColorBrewer")

col_fun = colorRamp2(c(0, 10,20), c("#081D58", "#41B6C4",  "#FFFF33" ))
color_info = c(rep("#000000", 2), rep("#1B9E77", 2), rep("#7570B3",4),
                rep("#000000",1), rep("#D95F02",2),rep("#000000",4),
                rep("#E6AB02",2),rep("#000000",7))

#
abundance_matrix <- as.matrix(abundance_order)

# row label 

row_labels_fdr_CVSP <-  as.vector(do.call(rbind,lapply(rownames(abundance_matrix), paste,"*",sep="")))

row_labels_fdr_CVSP[1] <- rownames(abundance_matrix)[1]
row_labels_fdr_CVSP[14] <- rownames(abundance_matrix)[14]
row_labels_fdr_CVSP[22] <- rownames(abundance_matrix)[22]
row_labels_fdr_CVSP[23] <- rownames(abundance_matrix)[23]
row_labels_fdr_CVSP[24] <- rownames(abundance_matrix)[24]

row_ha = rowAnnotation("Healthy vs ICU-" = group$`CVSM Pvalue sig`, 
                           "Healthy vs ICU+" = group$`CVSP Pvalue sig`,
                           "ICU- vs ICU+" = group$`MVSP Pvalue sig`,
                           col = list("Healthy vs ICU-" = c("Yes" = "black","No" ="grey"),
                                      "Healthy vs ICU+" = c("Yes" = "black","No" ="grey"),
                                      "ICU- vs ICU+" = c("Yes" = "black","No" ="grey")),
                       show_legend = FALSE,
                       simple_anno_size = unit(3, "mm"),
                       annotation_name_gp = gpar(fontsize=8.5,fontface = "bold"))

p <- Heatmap(abundance_matrix,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_labels = row_labels_fdr_CVSP,
        name = "Abundance",
        # mark colors 
        col = col_fun,
        # row_names_gp = gpar(col = c(rep(col[1], 13), rep(col[2], 8), rep(col[3],3))),
        row_names_gp = gpar(col = color_info,fontface = "italic"),
        #row_names_gp = gpar(col =  c(rep("red", 12), rep("blue", 12))),
        # split based on column 
        column_split = factor(c(rep("Healthy",11),rep("ICU-",14),rep("ICU+",49))),
        column_title = c("", "", ""),
        left_annotation = row_ha,
        #column_title_gp = gpar(fill = c("red", "blue", "green")),
        border = TRUE,
        top_annotation = HeatmapAnnotation(Group = c(rep("Healthy",11),rep("ICU-",14),rep("ICU+",49)),  col = list(Group = c("Healthy" = "#7FCDBB", "ICU-" = "#1D91C0", "ICU+" = "#081D58"))))#column_km = 3


lgd_sig = Legend(labels = c("Yes","No"), title = "Significance", legend_gp = gpar(fill = c("black","grey")),ncol = 1)

draw(lgd_sig)

pdf("result/significantspecies.pdf",width = 14)        
draw(p,heatmap_legend_side = "left", annotation_legend_side = "left", 
     annotation_legend_list = list(lgd_sig), merge_legend = TRUE)

grid.text(bquote(underline("Genus")), x=0.005, y=0.3,  gp=gpar(fontsize=10, col=c("#000000")), just = "left")

grid.text("Candida", x=0.005, y=0.275,  gp=gpar(fontsize=10, col=c("#7570B3")), just = "left")
grid.text("Aspergillus", x=0.005, y=0.25,  gp=gpar(fontsize=10, col=c("#1B9E77")), just = "left")

grid.text("Debaryomyces", x=0.005, y=0.225,  gp=gpar(fontsize=10, col=c("#D95F02")), just = "left")
grid.text("Penicillium", x=0.005, y=0.2,  gp=gpar(fontsize=10, col=c("#E6AB02")), just = "left")

grid.text("Others", x=0.005, y=0.175,  gp=gpar(fontsize=10, col=c("#000000")), just = "left")


grid.roundrect(x = 0.003, y = 0.235, width= 0.07, height=0.17,  name="box", just = "left")
dev.off()
###########

#unique(str_split(rownames(abundance_matrix)," ",simplify = T)[,1]
#Alternaria        Arxiella     Aspergillus         Candida      Chaetomium 
#1               1               2               4               1 
#Debaryomyces      Dipodascus  Geminibasidium      Gibberella Holtermanniella 
#2               1               1               1               1 
#Penicillium          Pichia   Saccharomyces       Saitozyma   Schizothecium 
#2               1               1               1               1 
#Serendipita    Trichosporon   Vishniacozyma 
#1               1               1 
#(candida;Aspergillus;Debaryomyces;Penicillium:i think these four are enough)




