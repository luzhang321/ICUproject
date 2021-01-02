# for ICU ITS 
# order : alpha diversity analysis; beta diversity analysis; top phylum analysis; 
# name : cuiqin 
# aim :  diversity & top phylum analysis 

setwd("/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS")
setwd("/media/lu/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS")
library(metagenomeSeq)
library(vegan)
library(OTUtable)
library(data.table)
library(ggpubr)
library(tidyverse)

# normalization --------------------------------------------------------------------

table_species_DADA2 <- fread("table_species_DADA2_human.csv", header = T)
#load("dada_output.RData")
table_species_DADA2 <- table_species_DADA2 %>%
  column_to_rownames(.,var="species") %>%
  .[,c(2:75)]

dada_norm <- table_species_DADA2 %>%   
  metagenomeSeq::newMRexperiment (.) %>%       ## creating MR object                      
  metagenomeSeq::cumNorm(., p = metagenomeSeq::cumNormStatFast(.)) %>% ## normalization  
  metagenomeSeq::MRcounts (., norm = TRUE, log = TRUE) %>%
  as_tibble(.) %>%
  add_column(otuid = rownames(table_species_DADA2)) %>% 
  column_to_rownames(.,var = "otuid") %>%
  filter_taxa(., abundance = 0, persistence = 10)


dada_diversity_input <- t(dada_norm)  %>%   
  data.frame()

write.table(dada_norm, "result/dada_species_10per.txt",row.names = T)  




# group information ---------------------------------------------------------------------
# 2.1 prepare group information
#ICU_pheno <- fread(file = "/Volumes/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt",header = T)
ICU_pheno <- fread(file = "/media/lu/Lucy/documents/201910_ICU_ITS/02_Normalization/LOCALTEST/ICU_pheno_anno.txt",header = T)
unique(ICU_pheno$group)
#[1] "Pip"      "Mero"     "ICU-"     "Control "
ICU_pheno$group <- factor(ICU_pheno$group, levels = c("Control", "ICU-", "ICU+Mero","ICU+Pip"))
ICU_pheno$Bgroup <- factor(ICU_pheno$Bgroup, levels = c("Control", "ICU-", "ICU+"))
ICU_pheno <- ICU_pheno[-c(74),] # remove AT081,not sequenced in ITS 


# alpha diversity ---------------------------------------------------------------------

# alpha diversity - shannon index 
ICU_norm_shannon <- diversity(dada_diversity_input, "shannon")
# alpha diversity - simpson index
ICU_norm_simpson <- diversity(dada_diversity_input, "simpson")
# alpha diversity - chao1 index 
ICU_norm_chao1 <- data.frame(sample = rownames(dada_diversity_input))
for (i in (1:nrow(dada_diversity_input))){
  ICU_norm_chao1[i,2] <- fossil::chao1(dada_diversity_input[i,])
}
colnames(ICU_norm_chao1)[2] <- "chao1"

# build a whole data.frame 
ICU_data_df <- data.frame(shannon = ICU_norm_shannon, simpson = ICU_norm_simpson, chao1 = ICU_norm_chao1$chao1) %>%
  data.frame(sample = rownames(.)) %>%
  merge(., ICU_pheno, ID = sample, sort = FALSE)


# plot 
kruskal_plot <- function(df, groupinfo, index, colorinfo, comparisons){
  ggboxplot(df, x = groupinfo, y = index, width = 0.6, color = groupinfo, add = "jitter", size = 1.6, add.params = list(size = 1.8), xlab = "Group", palette = colorinfo) +
    border() + # add border line 
    stat_compare_means(method = "kruskal.test", label.x.npc = "left", size = 5.7, label.y = 0.6) + 
    stat_compare_means(comparisons = comparisons) + 
    font("xlab", size = 12, face = "bold", color = "#2D2D2D") + 
    font("xy.text", size = 10, color = "#2D2D2D", face = "bold") +
    font("ylab", size = 12, face = "bold", color = "#2D2D2D") + 
    scale_fill_manual(values = alpha(colorinfo, .6)) + 
    main_theme + 
    #stat_compare_means(label = "p.signif", ref.group = "Control", size = 6) + # Pairwise comparison against Control
    rremove("legend") #remove the legend 
}

main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=10),
                   axis.title=element_text(size=10,face="bold"),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=8),
                   text=element_text(family="sans", size=8),
                   plot.title = element_text(hjust = 0.5,face="bold"))

ICU_data_df$group <- str_replace_all(ICU_data_df$group,"ICU\\+Pip","Pip") 
ICU_data_df$group <- str_replace_all(ICU_data_df$group, "ICU\\+Mero", "Mero")
ICU_data_df$group <- str_replace_all(ICU_data_df$group, "Control", "Healthy")
ICU_data_df$Bgroup <- str_replace_all(ICU_data_df$Bgroup, "Control", "Healthy")

ICU_data_df$group <- factor(ICU_data_df$group, levels = c("Healthy","ICU-","Mero","Pip"))
ICU_data_df$Bgroup <- factor(ICU_data_df$Bgroup, levels = c("Healthy","ICU-","ICU+"))

color = c("#f6921d", "#ed3224","#35454b","purple")
my_comparisons <- list(c("ICU-", "Healthy"), c("Mero","Healthy"), c("Pip", "Healthy"), c("Mero","ICU-"), c("Pip", "ICU-"),  c("Pip","Mero"))
k1_shannon <- kruskal_plot(ICU_data_df, "group", "shannon", color,my_comparisons)
k1_shannon

k1_simpson <- kruskal_plot(ICU_data_df, "group", "simpson", color,my_comparisons)
k1_chao1 <- kruskal_plot(ICU_data_df, "group", "chao1", color,my_comparisons)

color = c("#f6921d", "#ed3224", "#3b57a6")
my_comparisons <- list(c("ICU+", "ICU-"), c("ICU-", "Healthy"), c("ICU+", "Healthy") )
k2_shannon <- kruskal_plot(ICU_data_df, "Bgroup", "shannon", color,my_comparisons)
k2_simpson <- kruskal_plot(ICU_data_df, "Bgroup", "simpson", color,my_comparisons)
k2_chao1 <- kruskal_plot(ICU_data_df, "Bgroup", "chao1", color,my_comparisons)

p_alpha_diversity <- ggarrange(k1_shannon, k1_simpson, k1_chao1, k2_shannon, k2_simpson, k2_chao1, labels = c("A.","B.","C.","D.","E.","F."),  ncol = 3, nrow = 2)

# increase font size 10.31 
k1_shannon <- k1_shannon + font("xylab", size = 16) +
  font("xy.text", size = 16)
k1_simpson <- k1_simpson + font("xylab", size = 16) +
  font("xy.text", size = 16)
k1_chao1 <- k1_chao1 + font("xylab", size = 16) +
  font("xy.text", size = 16)

k2_shannon <- k2_shannon + font("xylab", size = 16) +
  font("xy.text", size = 16)
k2_simpson <- k2_simpson + font("xylab", size = 16) +
  font("xy.text", size = 16)
k2_chao1 <- k2_chao1 + font("xylab", size = 16) +
  font("xy.text", size = 16)

p_alpha_diversity <- ggarrange(k1_shannon, k1_simpson, k1_chao1, k2_shannon, k2_simpson, k2_chao1,ncol = 3, nrow = 2)

ggsave(p_alpha_diversity, filename = "result/S1.alpha_diversity_COLOR.png", width = 12,height = 10)
save.image("result/202006_latest.RData")
pdf("result/S1.alpha_diversity_COLOR.pdf", width = 13,height = 10)

p_alpha_diversity

dev.off()



# plot for beta diversity 
set.seed(2)
ICU_sel_bray_curtis <- vegdist(dada_diversity_input, method = "bray", diag = T, upper = T)
pcoa <- cmdscale(ICU_sel_bray_curtis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
# add sample infomation 
points <- points[sort(rownames(points)),]
points <- cbind(points, ICU_pheno)

# plot function 

library(RColorBrewer)
# draw the plot info =============================================================================================
# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
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

ICU_PcoA_plots <- function(points_df, matrix_t, groupinfo1, groupinfo2, title_name){
  ICU_pheno_sel <- subset(points_df, Bgroup %in% c(groupinfo1, groupinfo2))
  ICU_mat_t_sel <- subset(matrix_t, rownames(matrix_t) %in% as.character(ICU_pheno_sel$sample))
  ICU_mat_t_sel <- ICU_mat_t_sel[sort(rownames(ICU_mat_t_sel)),]
  set.seed(4)
  adonis_result <- adonis(ICU_mat_t_sel ~ Bgroup, data = ICU_pheno_sel, distance = "bray")
  pvalue <- adonis_result$aov.tab$`Pr(>F)`[1]
  Fvalue <- round(adonis_result$aov.tab$F.Model[1], 3)
  print(adonis_result)
  p <- ggplot(ICU_pheno_sel, aes(x=x, y=y, color=Bgroup))
  p <- p + geom_point(alpha=.7, size=1.5) +
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
         title = title_name) + main_theme +
    #annotate(geom = "text", label = paste("p-value",pvalue, sep = ":"),x=(max(ICU_pheno_sel$x)), y = min(ICU_pheno_sel$y),hjust=1,color = "black", size = 3) 
    annotate(geom = "text", label = paste("F:",Fvalue, "; pvalue:",pvalue, sep = " "),  x=max(points_df$x), y = min(points_df$y),hjust=1,  color = "black", size = 3)
  
  q <- p #scale_color_brewer(palette = "Set1")
  
  return(q)
}


pcoa1 <- ICU_PcoA_plots(points, dada_diversity_input, "Control", "ICU+","Healthy vs ICU+") + scale_color_manual(values = c("#f6921d", "#3b57a6"),labels = c("Healthy","ICU+"),name = "Group")

pcoa1

pcoa2 <- ICU_PcoA_plots(points, dada_diversity_input, "Control","ICU-","Healthy vs ICU-") +  scale_color_manual(values = c("#f6921d", "#ed3224"),labels = c("Healthy","ICU-"),name = "Group")

pcoa3 <- ICU_PcoA_plots(points, dada_diversity_input, "ICU-", "ICU+","ICU- vs ICU+") + scale_color_manual(values = c("#ed3224", "#3b57a6"),labels = c("ICU-","ICU+"),name = "Group")

ggarrange(pcoa1,pcoa2,pcoa3,ncol = 3,nrow=2)
 

#phylum 


# top3 phylum 
species_phylum <- read.csv("species_matched_phylum.csv",sep = ";") %>%
  .[,2:3] %>%
  filter(species !="NANA") %>%
  unique(.)


phylum_species_input <- dada_norm %>%
  mutate(species = rownames(dada_norm))

phylum_species_combine  <- merge(species_phylum, phylum_species_input, id = "species") %>% 
  select(.,-species) %>%
  group_by(., phylum) %>%
  summarise_all(., sum) %>%
  column_to_rownames(.,var="phylum")

phy_sum <- data.frame(Physum = rowSums(phylum_species_combine)) %>% 
  add_column(.,Phy = rownames(.)) %>%
  .[order(.$Physum, decreasing = TRUE),]

#Physum                      Phy
#p__Ascomycota            12687.05684            p__Ascomycota
#p__Basidiomycota          1457.17308         p__Basidiomycota
#p__Mortierellomycota       509.93884     p__Mortierellomycota
#p__Mucoromycota            107.09407          p__Mucoromycota
#p__Neocallimastigomycota    54.06661 p__Neocallimastigomycota

# calculate percentage 

write_delim(phy_sum, "phylum_species_combine.txt",delim = "\t")

summm <- sum(phy_sum$Physum) #14815.33
#p__Ascomycota

(12687.05684/summm)*100 # [1] 85.63466

# p__Basidiomycota
(1457.17308/summm)*100 # 9.835577

# p__Mortierellomycota

(509.93884/summm)*100 #[1] 3.441968

# normalization 
zero_to_one_normalization <- function(vector){
  max_num <- max(vector)
  min_num <- min(vector)
  nor_vec <- c()
  for (i in 1:length(vector)){
    nor_vec[i] <- (vector[i] - min_num)/(max_num - min_num)
  }
  return(nor_vec)
}


### for Ascomycota 
Ascomycota_abundance <- data.frame(phylum_species_combine[1,]) %>%
  .[,sort(colnames(.))] %>%
  t(.) %>% 
  as.data.frame(.) %>%
  rownames_to_column(.,var="sample") %>%
  add_column(.,nor=zero_to_one_normalization(.[,2])) 


### For Basidiomycota

Basidiomycota_abundance <- data.frame(phylum_species_combine[2,]) %>%
  .[,sort(colnames(.))] %>%
  t(.) %>% 
  as.data.frame(.) %>%
  rownames_to_column(.,var="sample") %>%
  add_column(.,nor=zero_to_one_normalization(.[,2])) 


# FOR  Mortierellomycota===============
Mortierellomycota_abundance <- data.frame(phylum_species_combine[3,]) %>%
  .[,sort(colnames(.))] %>%
  t(.) %>% 
  as.data.frame(.) %>%
  rownames_to_column(.,var="sample") %>%
  add_column(.,nor=zero_to_one_normalization(.[,2])) 

# for PcoA =============================
# pcoa this pcoa is built by previous with all selective species 
# 
pcoa_plot_rabinow <- function(data_pcoa, pheno_info, abundance_info, x, y, abundance, phylum){
  data_points_plot <- as.data.frame(data_pcoa$points)
  colnames(data_points_plot) <- c("x", "y", "z")
  eig <- data_pcoa$eig
  data_points_plot <- data_points_plot[sort(rownames(data_points_plot)),] %>%
    cbind(., pheno_info, abundance_info[,2:3])
  colnames(data_points_plot)[8] <- "Abundance"
  p <- ggplot(data_points_plot, aes(x=x, y=y, color=Abundance, shape = Bgroup))
  p <- p + geom_point(alpha=.7, size=1.5) +
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
         title=phylum) + main_theme +
    scale_shape_discrete(name = "Group",labels = c("Healthy","ICU-","ICU+")) 
  pcoaplot <- p + scale_color_gradientn(colours=rainbow(7))
  return(pcoaplot)
}

q1 <- pcoa_plot_rabinow(pcoa, ICU_pheno, Ascomycota_abundance, x, y, abundance, "Ascomycota") + theme(legend.key.size = unit(0.3, "cm"))
q2 <- pcoa_plot_rabinow(pcoa, ICU_pheno, Basidiomycota_abundance, x, y, abundance, "Basidiomycota") + theme(legend.key.size = unit(0.3, "cm"))
q3 <- pcoa_plot_rabinow(pcoa, ICU_pheno, Mortierellomycota_abundance, x, y, abundance, "Mortierellomycota") + theme(legend.key.size = unit(0.3, "cm"))
q1
ggarrange(q1,q2,q3)


p_pcoa_phylum_beta <- ggarrange(pcoa1, pcoa2, pcoa3, q1, q2, q3, nrow = 2, ncol = 3)
write(pcoa$points, "result/PCoA.result.txt")
ggsave("result/S3.PCoA.Phylum_Beta.png",p_pcoa_phylum_beta,width = 10,height = 5)

pdf("result/S3.PCoA.Phylum_Beta.pdf",width = 10,height = 5)
p_pcoa_phylum_beta
dev.off()

save.image("result/202007.ITS.phylum.diversity.RData")
load("result/202007.ITS.phylum.diversity.RData")

