# Metabolites
# metabolites ~ sig.diff bacteria species correlation 
# name : cuiqin 
# sig.diff : higher in healthy vs icu+(not in icu-); higher in icu- vs icu+
# cor : spearman correlation

library(data.table)
library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(readxl)
library("RColorBrewer")
library(circlize)
library(pacman)
p_load(ppcor)
setwd("/media/lu/Lucy/documents/201910_ICU_ITS/08_metabolites_enrichmentspecies/")


# read file - species list 
# ----------------------------------------------------------------
species_list <- read_excel("/media/lu/Lucy/documents/201910_ICU_ITS/08_metabolites_enrichmentspecies/data/wilcox_diff_species_group.xlsx")

colnames(species_list) <- c("Control_ICUm","Higher1","Control_ICUp","Higher2","ICUm_ICUp","Higher3")

Control_ICUm_higher <- filter(species_list, Higher1 %in% "Control") %>%
  .$Control_ICUm

Control_ICUp_higher <- filter(species_list, Higher2 %in% "Control") %>%
  .$Control_ICUp

Control_both_higher <- setdiff(Control_ICUp_higher,Control_ICUm_higher) # the species only in Control vs ICU+ : 31 of them, one is removed later, cause it is only in one sample, and it's impossible to do correlation
# ignore the name of vector

ICUm_higher <- filter(species_list, Higher3 %in% "ICU-") %>% # 8 
  .$ICUm_ICUp


intersect(Control_both_higher, ICUm_higher) # no intersection 

sig_list <- union(Control_both_higher, ICUm_higher) %>%
  str_replace_all(.," ","_") # in total 39 sig list 
  


#write the sig list
write_delim(data.table(sig_list = sig_list),"result/sig_species_correlation/change_species/change_species_2_healthy_only_higher_ICUp/sig_species_plot.txt")


# # read file - species abundance table 
# ----------------------------------------------------------------

abundance_table <- fread("data/ICU_merged_table_species.txt")

# read file - scfa 
# ----------------------------------------------------------------

scfa <- readxl::read_excel("data/SCFA_relative_value.xlsx")

# read file - ba 
# ----------------------------------------------------------------

ba <- readxl::read_excel("data/BA_quantified_value.xlsx")

# after sum up unknown bile acids in bile acids part, then it will end up with those like things -> unknown part, then all of them should be removed 
# only 7 bile acids are left 
# ----------------------------------------------------------------

sig_ba <- c("12-Ketolithocholic_acid","Deoxycholic_acid","Glycolithocholic_acid","Hyodeoxycholic_acid","Isolithocholic_acid","Lithocholic_acid","Ursodeoxycholic_acid") %>%
  str_replace_all(.,"_"," ")



sig_scfa <- c("Acetic_acid","Propionic_acid","Butyric_acid","Valeric_acid","Caproic_acid","Heptanoic_acid") %>%
  str_replace_all(.,"_"," ")  



# extract significant ones - scfa 
# ----------------------------------------------------------------

head(scfa)

scfa_sel <- scfa[,sig_scfa] %>%
  mutate(sample = scfa$sample)


# extract significant ones - bile acid
# ----------------------------------------------------------------

head(ba)

ba_sel <- ba[,sig_ba] %>%
  mutate(sample = ba$`Customer ID`)


# extract significant ones - species
# ----------------------------------------------------------------

abundance_table_sel <- abundance_table %>%
  column_to_rownames(.,var = "ID") %>%
  t(.) %>%
  as.data.frame() %>%
  rownames_to_column(.,var = "sample")


abundance_table_sel <- abundance_table_sel[,c(sig_list,"sample")] # extract the significant species list 


abundance_table_sel$sample <- str_remove_all(abundance_table_sel$sample,"_qc_noHuman_1_profile") # remove the _qc_noHuman_1_profile present in the sample name 
# only select the samples present in metabolites table 
abundance_table_sel_matrix <- abundance_table_sel %>%
  filter(., sample %in% scfa_sel$sample) %>%
  column_to_rownames(.,var = "sample") %>% # remove one column to rowname 
  as.matrix()


# make input for correlation - as.matrix 
# ----------------------------------------------------------------
scfa_sel_matrix <- scfa_sel %>%
  column_to_rownames(.,var = "sample") %>%
  as.matrix()

ba_sel_matrix <- ba_sel %>%
  column_to_rownames(.,var = "sample") %>%
  as.matrix()


# group information 

# group <- fread("/media/lu/Lucy/documents/201910_ICU_ITS/11_Grid/whole_dataset/Group_table.txt")
# group_sample <- data.table(sample = rownames(scfa_sel_matrix)) %>%
#   merge(., group, by = "sample", all.x = TRUE, sort = FALSE) %>%
#   mutate(Bgroup=replace(Bgroup,Bgroup == "ICU-", 1)) %>%
#   mutate(Bgroup=replace(Bgroup,Bgroup == "ICU+", 1)) %>%
#   mutate(Bgroup=replace(Bgroup,Bgroup == "Healthy", 0))


# correlation function (from mano,has been looked into details,good function)
# ----------------------------------------------------------------

Big <- abundance_table_sel_matrix
Small <- scfa_sel_matrix

# correlation function : make samples as the row names 
same_samples_Correlation <- function(Big,Small) {
  
  ####Corrrelation_Function
  Data_1 <- Big
  Data_2 <- Small
  #Data_1 = abundance_table_sel_matrix
  #Data_2 = scfa_sel_matrix
  
  Data_1 <- Data_1 %>% data.matrix()
  Data_2 <- Data_2 %>% data.matrix()
  
  Correlation_Matrix = matrix(nrow = ncol(Data_2),ncol = 1)
  #Correlation_Vector = vector()
  
  PValue_matrix = matrix(nrow = ncol(Data_2),ncol = 1)
  #Pvalue_vector = vector()
  
  for (i in 1:ncol(Data_1)) {
    Correlation_Vector = vector()
    Pvalue_vector = vector()
    for (j in 1:ncol(Data_2)) {
      # which method to use ? pcor.test from pcor package 
      Correlation = cor.test(x = Data_1[,i], y = Data_2[,j],method = "spearman")
      
      Current_line = Correlation$estimate
      Current_pvalue = Correlation$p.value
      Correlation_Vector = rbind(Correlation_Vector,Current_line)
      Pvalue_vector = rbind(Pvalue_vector,Current_pvalue)
      
    }
    Correlation_Matrix = cbind(Correlation_Matrix,Correlation_Vector)
    PValue_matrix = cbind(PValue_matrix,Pvalue_vector)
  }
  
  Correlation_Matrix = Correlation_Matrix[,-1] # first column is NA 
  PValue_matrix = PValue_matrix[,-1] # first column is NA 
  
  rownames(Correlation_Matrix) = colnames(Data_2)
  rownames(PValue_matrix) = colnames(Data_2)
  
  colnames(Correlation_Matrix) = colnames(Data_1)
  colnames(PValue_matrix) = colnames(Data_1)
  
  
  
  QValue_matrix = matrix(p.adjust(as.vector(as.matrix(PValue_matrix)), method='fdr'),ncol=ncol(PValue_matrix))  
  rownames(QValue_matrix) = rownames(PValue_matrix)
  colnames(QValue_matrix) = colnames(PValue_matrix)
  
  Qvalue_bh_matrix = matrix(p.adjust(as.vector(as.matrix(PValue_matrix)), method='BH'),ncol=ncol(PValue_matrix))  
  rownames(Qvalue_bh_matrix) = rownames(PValue_matrix)
  colnames(Qvalue_bh_matrix) = colnames(PValue_matrix)
  
  
  list = list("cor" = Correlation_Matrix,"pvalue" = PValue_matrix,"qvalue"= QValue_matrix,"qvalue-BH" = Qvalue_bh_matrix)
  return(list)
  
}


# calculate spearman correlation 
Summary_Species_scfa <- same_samples_Correlation(abundance_table_sel_matrix, scfa_sel_matrix)


Summary_Species_ba <- same_samples_Correlation(abundance_table_sel_matrix, ba_sel_matrix)



# correlation combine 
# ----------------------------------------------------------------

correlation <- rbind(Summary_Species_scfa$cor,Summary_Species_ba$cor) %>%
  as.matrix()
pvalue <- rbind(Summary_Species_scfa$pvalue,Summary_Species_ba$pvalue) %>%
  as.matrix()
qvalue <- rbind(Summary_Species_scfa$qvalue,Summary_Species_ba$qvalue) %>%
  as.matrix()

qvalue_bh <- rbind(Summary_Species_scfa$`qvalue-BH`,Summary_Species_ba$`qvalue-BH`) %>%
  as.matrix()
# write the correlation to file 
# ----------------------------------------------------------------

correlation_write <- correlation %>%
  as.data.frame() %>%
  rownames_to_column(.)

write_delim(correlation_write,"result/sig_species_correlation/change_species/change_species_2_healthy_only_higher_ICUp/correlation.csv",delim = "\t")

qvalue_write <- qvalue %>%
  as.data.frame() %>%
  rownames_to_column(.)

write_delim(qvalue_write,"result/sig_species_correlation/change_species/change_species_2_healthy_only_higher_ICUp/correlation-qvalue.csv",delim = "\t")

pvalue_write <- pvalue %>%
  as.data.frame() %>%
  rownames_to_column(.)

write_delim(pvalue_write,"result/sig_species_correlation/change_species/change_species_2_healthy_only_higher_ICUp/correlation-pvalue.csv",delim = "\t")

# read file of GRiD 
# ----------------------------------------------------------------------------------------

GRiD <- fread("/media/lu/Lucy/documents/201910_ICU_ITS/11_Grid/whole_dataset/merged_table.txt") %>%
  filter(Genome %in% sig_list) %>%
  merge(data.table(Genome = sig_list), ., by = "Genome", all.x = TRUE, sort = F)

GRiD_ann <- GRiD %>%
  column_to_rownames(., var = "Genome") %>%
  t(.) 

GRiD_ann <- GRiD_ann[,-10]
GRiD_ann_bar <- HeatmapAnnotation(GRiD = anno_boxplot(GRiD_ann, height = unit(1.5, "cm"), 
                                                      size = unit(1, "mm")),
                                  show_annotation_name = TRUE)



# sig species mark 

ICUm_higher <- filter(species_list, Higher3 %in% "ICU-") %>% # 8 
  .$ICUm_ICUp

sig_mark <- ifelse(colnames(correlation) %in% ICUm_higher, "ICU-","Healthy")
# ifelse(colnames(correlation) %in% "Bacteroides_caccae")
table(sig_mark) # 31(Healthy vs ICU+ (but not in ICU-)) 8(ICU- vs ICU+) 


# fva metabolites model ( because of changeing species, some might not be tested )
# -----------------------------------------------------------------------------

# merge it with metadata 

species_name <- data.table(species = colnames(correlation)) 
species_name$second <- str_split(species_name$species, "_", simplify = T)[,2]
table(species_name$second)# the ones with replicated are : bacterium :5;sp:3; unclassfied 5.

# the one with value : secretion value  

fva_matrix_scfa_icuEc_atp_wrap_to_merge_ori <- fread("result/sig_species_correlation/change_species/change_species_2_healthy_only_higher_ICUp/fva_matrix_scfa_icuEc_atp_wrap_V4.csv")


# now for the species has mulitple matched names, i will choose the highest one, for the ones with only one matched name,
# i will use that one. 
fva_matrix_scfa_icuEc_atp_wrap_to_merge <- fva_matrix_scfa_icuEc_atp_wrap_to_merge_ori %>%
  add_column(second = str_split(.$V1, "_", simplify = T)[,2]) %>%
  left_join(species_name, . , by = "second") %>%
  filter(!(species == "Bacteroides_sp_9_1_42FAA" & !V1 == "Bacteroides_sp_9_1_42FAA")) %>%
  filter(!(species == "Clostridium_sp_ATCC_BAA_442" & !V1 == "Clostridium_sp_L2_50")) %>%
  filter(!(species == "Ruminococcus_sp_JC304" & !V1 == "Ruminococcus_sp_5_1_39BFAA")) %>% 
  mutate(V1 = replace(V1, species == "Bilophila_unclassified", "Bilophila_unclassified")) %>%
  mutate(V1 = replace(V1, species == "Coprobacter_fastidiosus", "Coprobacter_fastidiosus")) %>%
  mutate(V1 = replace(V1, species == "Megasphaera_unclassified", "Megasphaera_unclassified")) %>%
  mutate(V1 = replace(V1, species == "Oscillibacter_unclassified", "Oscillibacter_unclassified")) %>%
  mutate(V1 = replace(V1, species == "Subdoligranulum_unclassified", "Subdoligranulum_unclassified")) %>%
  mutate(V1 = replace(V1, species == "Veillonella_unclassified", "Veillonella_unclassified")) %>%
  #!!!!! these ones will be summarized later  Lachnospiraceae_bacterium_1_1_57FAA; Lachnospiraceae_bacterium_2_1_58FAA
  filter(!(species == "Lachnospiraceae_bacterium_1_1_57FAA" & !V1 == "Lachnospiraceae_bacterium_sp_8_1_57FAA")) %>% 
  filter(!(species == "Lachnospiraceae_bacterium_2_1_58FAA" & !V1 == "Lachnospiraceae_bacterium_sp_8_1_57FAA")) %>% 
  # filter(!(species == "Lachnospiraceae_bacterium_1_1_57FAA" & !V1 == "Lachnospiraceae_bacterium_1_1_57FAA")) %>%
  #!!!!! these ones will be summarized later Erysipelotrichaceae_bacterium_21_3;Erysipelotrichaceae_bacterium_2_2_44A;
  # Erysipelotrichaceae_bacterium_6_1_45
  filter(!(species == "Erysipelotrichaceae_bacterium_21_3" & !V1 == "Erysipelotrichaceae_bacterium_5_2_54FAA")) %>% 
  filter(!(species == "Erysipelotrichaceae_bacterium_2_2_44A" & !V1 == "Erysipelotrichaceae_bacterium_5_2_54FAA")) %>% 
  filter(!(species == "Erysipelotrichaceae_bacterium_6_1_45" & !V1 == "Erysipelotrichaceae_bacterium_5_2_54FAA")) %>% 
  dplyr::select(., -V1, -second, -obj) %>%
  pivot_longer(., -species, names_to = "SCFA", values_to = "Secrection") %>%
  pivot_wider(., names_from = "species", values_from = "Secrection") %>%
  column_to_rownames(., var = "SCFA") %>%
  as.matrix()
  

fva_matrix_scfa_icuEc_atp_wrap_to_merge <- fva_matrix_scfa_icuEc_atp_wrap_to_merge[,-10] # remove one species because this species doens't have correlation value 
View(fva_matrix_scfa_icuEc_atp_wrap_to_merge)
colnames(fva_matrix_scfa_icuEc_atp_wrap_to_merge) 



# the one without value : only secretion or not 


fva_matrix_scfa_icuEc_atp_bin_wrap_to_merge <- fva_matrix_scfa_icuEc_atp_wrap_to_merge 

fva_matrix_scfa_icuEc_atp_bin_wrap_to_merge[fva_matrix_scfa_icuEc_atp_bin_wrap_to_merge>0] <- 1 

View(fva_matrix_scfa_icuEc_atp_bin_wrap_to_merge)

# compplexheatmap small part1 

#col_metabolites_model = colorRamp2(c(0,0.25,0.5,1), c("#4092b9", "#7abca2", "#b8e195", "#f6fbb5"))
col_metabolites_model = colorRamp2(c(0,0.2,0.5,1), c("#428ab7", "#7abca2", "#b8e195", "#f3faad"))
metabolites_model_combine <- ComplexHeatmap::Heatmap(fva_matrix_scfa_icuEc_atp_wrap_to_merge,
                        col = col_metabolites_model,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        show_heatmap_legend = FALSE,
                        na_col = "#979e99", # if some value is NA, will replaced by this color   #eff0eb
                        border = TRUE,#8. border 
                        row_dend_reorder = TRUE,
                        column_labels = str_replace_all(colnames(fva_matrix_scfa_icuEc_atp_bin_wrap_to_merge),"_"," "),
                        column_names_max_height=unit(8,"cm"),
                        column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                        row_names_gp = gpar(fontsize = 10),
                        left_annotation = rowAnnotation(Metabolites = c(rep("SCFA",3)), col = list(Metabolites = c("SCFA" = "#F7B22A")),border = TRUE, show_annotation_name = FALSE,simple_anno_size = unit(0.28, "cm"), annotation_legend_param = list(Metabolites = list(title = " ")), annotation_name_gp = gpar(fontsize=10),show_legend = F),# 4.left annotation  c("SCFA" = "#157CE0", "BA" = "#AF2912")
                        row_title = "Modeling",
                        row_title_gp = gpar(fontsize = 10)
                        # layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c) {
                        #   
                        #   ind_mat = restore_matrix(j, i, x, y)
                        #   ind = unique(c(ind_mat[1, 1]))
                        #   print(ind)
                        #   grid.rect(x[ind], y[ind], gp = gpar(lwd = 2, fill = "transparent"))
                        # }
                        )
draw(metabolites_model_combine)
# metabolites model yes or no 
#col_metabolites_model_yes_no = structure(c("#428ab7","#f37c54"), names = c("0", "1"))

col_metabolites_model_yes_no = structure(c("#428ab7","#f3faad"), names = c("0", "1"))

# col_metabolites_model = colorRamp2(c(0,0.2,0.5,1), c("#428ab7", "#7abca2", "#b8e195", "#f3faad"))
metabolites_model_combine_yes_no <- ComplexHeatmap::Heatmap(fva_matrix_scfa_icuEc_atp_bin_wrap_to_merge,
                                                     col = col_metabolites_model_yes_no,
                                                     cluster_columns = FALSE,
                                                     cluster_rows = FALSE,
                                                     show_heatmap_legend = FALSE,
                                                     na_col = "#979e99", # if some value is NA, will replaced by this color   #eff0eb
                                                     border = TRUE,#8. border 
                                                     row_dend_reorder = TRUE,
                                                     column_labels = str_replace_all(colnames(fva_matrix_scfa_icuEc_atp_bin_wrap_to_merge),"_"," "),
                                                     column_names_max_height=unit(8,"cm"),
                                                     column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                                                     row_names_gp = gpar(fontsize = 10),
                                                     left_annotation = rowAnnotation(Metabolites = c(rep("SCFA",3)), col = list(Metabolites = c("SCFA" = "#F7B22A")),border = TRUE, show_annotation_name = FALSE,simple_anno_size = unit(0.28, "cm"), annotation_legend_param = list(Metabolites = list(title = " ")), annotation_name_gp = gpar(fontsize=10),show_legend = F),# 4.left annotation  c("SCFA" = "#157CE0", "BA" = "#AF2912")
                                                     row_title = "Metabolic modeling",
                                                     row_title_gp = gpar(fontsize = 10)
                                                     # layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c) {
                                                     #   
                                                     #   ind_mat = restore_matrix(j, i, x, y)
                                                     #   ind = unique(c(ind_mat[1, 1]))
                                                     #   print(ind)
                                                     #   grid.rect(x[ind], y[ind], gp = gpar(lwd = 2, fill = "transparent"))
                                                     # }
)
draw(metabolites_model_combine_yes_no)



# complexheatmap 

# complex heatmap   -  correlation part 
# ----------------------------------------------------------------



# 1. test which color is situable 

col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c("#245b87", "#58b7e0","#FFFFFF","#C24933","#AF2912"))

#col_fun = colorRamp2(c(-1,0,1), c("#024C89","#FFFFFF","#BC0029"))



# 2. add annotation 
# : https://color.adobe.com/zh/create/color-wheel


# 3. add fixed qvalue 

# filter depends on correlation plot : below it is the version I tried with different filter, but finall
# we still use the first one ( only 0.05 no filter of correlation)

correlation
# fixed_qvalue <- qvalue # if the correlation <= 0.05, then i will not show it in the figure. 
# 
# 
# for (i in 1:nrow(correlation)){
#   for (j in 1:ncol(correlation)){
#     #print(testtt[i,j])
#     if (abs(correlation[i,j]) < 0.35){
#       fixed_qvalue[i,j] <- 1
#       #print(paste(i,j,sep = ","))
#       #print(correlation[i,j])
#     }
#   }
# }


# to highlight the species we are interested 

interested_matrix <- qvalue 
# 
# for (j in 1:ncol(correlation)){
#     for (i in 1:nrow(correlation)){
#       if (j %in% c(1,4,8,27,45,53)){
#         print(j)
#       }else{
#         interested_matrix[i,j] <- 1
#       }
#     }
#     
# }


# remove the Bacteroides_sp_9_1_42FAA : this species only exists in one sample, non-sense
correlation_new <- correlation[,-10]
qvalue_new <- qvalue[,-10] # real qvalue doesn't affect the result, when there are NAs or not 

# when I draw plot, I separately use the below 3 commands 


p <- Heatmap(correlation_new,
             name = "Spearman's rank correlation",
             col = col_fun,
             na_col = "#686868",
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #  if(qvalue[i, j] < 0.005)
             #    grid.circle(x = x, y = y, r = 0.02,gp = gpar(col = "white"))
             #},
             # label the circle of fdr 
             layer_fun = function(j, i, x, y, width, height, fill) {
               # here i is the row number from 1-6 and from 7-13
               # j is the the column number from 1-57 ( two list : the top part and the down part)
               # x coordinate of middle point of the cell which is measured in the viewport of the heatmap body.
               # y coordinate of middle point of the cell which is measured in the viewport of the heatmap body
               ind_mat = restore_matrix(j, i, x, y)  
               
               
               v = pindex(qvalue_new, i, j)
               l = v < 0.05
               ind_mat_vec = as.matrix(ind_mat)
               ind = ind_mat_vec[l]
               grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
               
               
               # v2 = pindex(interested_matrix, i, j)
               # l2 = v2 < 0.05
               # ind_mat_vec2 = as.matrix(ind_mat)
               # ind2 = ind_mat_vec2[l2]
               # grid.points(x[ind2], y[ind2], pch = 1, size = unit(3, "mm"), gp = gpar(col = "white",lwd=2))
               
             },
             cluster_rows = FALSE,
             cluster_columns = FALSE, 
             show_heatmap_legend = FALSE,
             #annotation_name_gp = gpar(fontsize=2)
             top_annotation = HeatmapAnnotation(Group = c(rep("Healthy",30),rep("ICU-",8)),  col = list(Group = c("Healthy" = "#3C93FA", "ICU-" = "#AF2912")),border = TRUE,annotation_name_gp = gpar(fontsize=10),simple_anno_size = unit(0.28, "cm"),annotation_legend_param = list(Group = list(title = " ")),show_legend = F), # 3.top annotation  c("Healthy" = "#F7B22A", "ICU" = "#D01B30")
             left_annotation = rowAnnotation(Metabolites = c(rep("SCFA",6),rep("BA",7)), col = list(Metabolites = c("SCFA" = "#F7B22A", "BA" = "#FA1846")),border = TRUE, show_annotation_name = FALSE,simple_anno_size = unit(0.28, "cm"), annotation_legend_param = list(Metabolites = list(title = " ")), annotation_name_gp = gpar(fontsize=10),show_legend = F),# 4.left annotation  c("SCFA" = "#157CE0", "BA" = "#AF2912")
             #column_split = factor(c(rep("Healthy",39),rep("ICU-Enriched",23))), #5. split the heatmap based on the title 
             column_title = "Enriched Species in comparison to ICU+ (n = 38)", #6.column title 
             column_gap = unit(0,"mm"),
             column_title_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 10, fontface = "italic"),
             #column_names_gp = gpar(fontsize = 10, fontface = "italic", col = sig_species_mark$Mark), # column color 
             row_split = factor(c(rep("SCFA",6),rep("BA",7))),
             row_title = "Metabolites (n = 13)",
             row_title_gp = gpar(fontsize = 10),
             row_gap = unit(3, "mm"), # 7.split gap setting 
             row_names_gp = gpar(fontsize = 10),
             border = TRUE,#8. border 
             heatmap_legend_param = list(legend_direction = "horizontal",title_position = "topcenter"),# lenged horizontal 
             row_dend_reorder = TRUE,
             #bottom_annotation = GRiD_ann_bar,
             column_labels = str_replace_all(colnames(correlation_new),"_"," "),
             column_names_max_height=unit(8,"cm"))


lgd_fdr = Legend(pch = 1, type = "points", labels = "FDR q < 0.05", legend_gp = gpar(col = "white",lwd=2),background = "black", title = " ")

#lgd_sig_species = Legend(pch = 16, type = "points", labels = c("ICU vs Healhty","ICU- vs ICU+","Both"), legend_gp = gpar(col = c("#B33C14","#BF801B","#6ECC93")),background = "white", title = "Significantly different species") # nrow=1 makes the lgend horizontal # title_position = "leftcenter", title_gp = gpar(fontface = "plain",fontsize =10)

lgd_group = Legend(labels = c("Healthy","ICU-"), title = "Enriched groups", legend_gp = gpar(fill = c("#3C93FA","#AF2912")),nrow = 1)

lgd_metabolites = Legend(labels = c("SCFA","Bile Acids"), title = "Metabolites", legend_gp = gpar(fill = c("#F7B22A","#FA1846")),nrow = 1)

lgd_abundance = Legend(col_fun = col_fun, title = "Spearman's rank correlation", border = "black",direction = "horizontal",legend_width = unit(4.8, "cm"))

lgd_secretion = Legend(col_fun = col_metabolites_model, title = "Secretion", border = "black",direction = "horizontal",legend_width = unit(4.8, "cm"))

lgd_no_model = Legend(labels = c("No model"), title = "   ", legend_gp = gpar(fill = c("#979e99")),nrow = 1)
#lgd_na = Legend(labels = c("NA"), title = "   ", legend_gp = gpar(fill = c("#686868")),nrow = 1)


#lgd_list = list(lgd_fdr,lgd_metabolites,lgd_group,lgd_abundance,lgd_secretion)

pd <- packLegend(lgd_fdr,lgd_metabolites,lgd_group,lgd_abundance,lgd_secretion, lgd_no_model,max_width = unit(29, "cm"), column_gap = unit(8, "mm"), row_gap = unit(1, "cm"),direction = "horizontal")
ComplexHeatmap::draw(p,annotation_legend_list = pd,heatmap_legend_side = "bottom", merge_legend = TRUE)





# combine - part3 


le = c(rep("Healthy",30),rep("ICU-",8))
names(le) <- str_replace_all(colnames(correlation_new),"_"," ")
input <- rbind(letters = le) 
rownames(input) <- ""
ht3 = Heatmap(input, 
              name = " ", 
              col = c("Healthy" = "white", "ICU-" = "white"),
              height = unit(0.000002, "cm"),
              width = unit(15, "mm"),
              show_heatmap_legend = FALSE,
              column_labels = str_replace_all(colnames(correlation_new),"_"," "),
              column_names_max_height=unit(8,"cm"),
              column_title_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10, fontface = "italic"))


#p_combine <-  p %v% metabolites_model_combine %v% GRiD_ann_bar %v% ht3



# color gold,turquoise1,purple for the box 
# for yes or no model 
# -------------------------------

metabolites_model_combine_yes_no

pdf("result/sig_species_correlation/change_species/change_species_2_healthy_only_higher_ICUp/combine_box_name_adjust_length_yes_no.pdf",width = 12,height = 9)

# combine plot 

p_combine_yes_no <-  p %v% metabolites_model_combine_yes_no %v% GRiD_ann_bar %v% ht3

lgd_secretion_yes_no = Legend(labels = c("Yes","No"), title = "Secretion", legend_gp = gpar(fill = c("#f3faad","#428ab7")),nrow = 1)

pd_yes_no <- packLegend(lgd_fdr,lgd_metabolites,lgd_group,lgd_abundance,lgd_secretion_yes_no, lgd_no_model,max_width = unit(29, "cm"), column_gap = unit(8, "mm"), row_gap = unit(1, "cm"),direction = "horizontal")


draw(p_combine_yes_no,annotation_legend_list = pd_yes_no,heatmap_legend_side = "bottom", merge_legend = TRUE)

decorate_annotation("GRiD", {
  grid.lines(x = unit(c(0,39),"native"),
             y = unit(c(1,1),"native"),
             gp = gpar(col = "red",
                       lty = "dashed")) # dotted 
})



dev.off()



save.image("result/sig_species_correlation/change_species/change_species_2_healthy_only_higher_ICUp/heatmap.RData")

load("result/sig_species_correlation/change_species/change_species_2_healthy_only_higher_ICUp/heatmap.RData")


#good for color picking https://coolors.co/3c93fa-af2912-5b6c5d-3b2c35-2a1f2d 
