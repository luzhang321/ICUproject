# Resistome Analysis/Profiling 

# difference among different groups : dabestr 

# moving from pvalue 
# this paper : https://www.researchgate.net/publication/333884529_Moving_beyond_P_values_data_analysis_with_estimation_graphics 

# dabestr method 

#===========set the parameter ============================================================

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Rscript dabestr_ARGs.R arg-RPKM-file group-file outdir outname")
} else if (length(args)==4) {
  print(paste("arg-RPKM-file:",args[1],",group-file,",args[2],",outdir:",args[3],",outname:",args[4],sep = ""))
} 


# library
#library(tidyverse)
library(dabestr)
library(data.table)
library(ggpubr)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
# read file 

#setwd("/Volumes/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/webserver_identity80_cov50/")

ARGs_RPKM <- read_delim(args[1], delim = "\t")
#ARGs_RPKM <- read_delim("data/merged_16S_subtypes_cov50_identity80_website_new.txt", delim = "\t")


colnames(ARGs_RPKM)[2:ncol(ARGs_RPKM)] <- str_split(colnames(ARGs_RPKM)[2:ncol(ARGs_RPKM)],pattern = ".deep", simplify = T)[,1]


ICU_pheno <- fread(args[2])
#ICU_pheno <- fread("/Volumes/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/tidy_directory/script/Resistome_composition/data/ICU_pheno_anno.txt")

ICU_pheno$Group <- ICU_pheno$Bgroup

ICU_pheno[Bgroup == "Control"]$Group <- "Healthy"

ICU_pheno$Group <- factor(ICU_pheno$Group, levels = c("Healthy","ICU-","ICU+"))

# filter 10% prevalance.

ARGs_RPKM <- ARGs_RPKM[rowSums(ARGs_RPKM[,2:ncol(ARGs_RPKM)] >0 )>=7,] # 


# long format 

ARGs_RPKM_long <- ARGs_RPKM %>%
  pivot_longer(-ID, names_to = "sample", values_to = "RPKM") %>% 
  merge(.,ICU_pheno, id = "sample")


# all ARGs name : 159 ARGs 
ARGs_name <- unique(ARGs_RPKM_long$ID)
length(ARGs_name)

# dabest 

a <- 0 
b <- 0

dabestr_list <- list()
dabestr_list_result <- list()
dabestr_list_result_sel <- list()
arg_sel <- vector()

p_h <- list() 
p_h_sel <- list() 
p_l <- list()
p_l_sel <- list()

# For ICU- vs ICU+ 

j <- 1 
for (i in 1:length(ARGs_name)){
  arg <- ARGs_name[i]
  unpaired <- filter(ARGs_RPKM_long, ID == arg) %>% 
    dabest(Group, RPKM, 
         idx = list(c("ICU-", "ICU+")),
         paired = FALSE
    )

  
  if (unpaired$result$difference > 0){
    a<- a+1 #134
    p_h[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"",sep = ""))
  }else if (unpaired$result$difference < 0){
    b <- b+1 #200
    p_l[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"",sep = "")) 
  }
  
  if ((unpaired$result$difference > 0 & unpaired$result$bca_ci_low >0) ||  (unpaired$result$difference < 0 & unpaired$result$bca_ci_high < 0)){
    arg_sel[j] <- arg
    dabestr_list_result_sel[[j]] <- unpaired$result[,c(1:13)]
    
    j <- j + 1 
  }
  if ((unpaired$result$difference > 0 & unpaired$result$bca_ci_low >0)){

    p_h_sel[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"",sep = ""), palette = c("#ed3224","#3b57a6"))
    #pdf(paste(arg,".pdf",sep = ""))
    pdf(paste(args[3],args[4],arg,"ICU+vsICU-.high.pdf",sep = ""))
    print(p_h_sel[[arg]])
    dev.off()
  }else if ((unpaired$result$difference < 0 & unpaired$result$bca_ci_high < 0)){
    p_l_sel[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"",sep = ""))
    pdf(paste(args[3],args[4],arg,"ICU+vsICU-.low.pdf",sep = ""))
    print(p_h_sel[[arg]])
    dev.off()
  }
  
  dabestr_list[[i]] <- unpaired
  dabestr_list_result[[i]] <- unpaired$result
  
  
}

dabestr_out <- do.call(rbind, dabestr_list_result) 
dabestr_out$arg <- ARGs_name

dabestr_out_sel <- do.call(rbind, dabestr_list_result_sel) 
dabestr_out_sel$arg <- arg_sel

a_sel <- do.call(rbind, dabestr_list_result_sel) 
a_sel$arg <- arg_sel

dim(dabestr_out_sel)


dabestr_out <- dabestr_out[order(dabestr_out$difference, decreasing = T),]
#dim(dabestr_out_sel)

write_delim(dabestr_out_sel, paste(args[3],args[4],"dabeset.ICU.output.txt",sep = ""), delim = "\t")
#write_delim(dabestr_out_sel, "dabestr//dabeset.ICU.output.txt", delim = "\t")


ICU_p_h_sel 
p_h_sel

print(a)
print(b)


p_h_sel_p <- ggarrange(plotlist = p_h_sel, ncol = 3, nrow = 2)
if (length(names(p_l_sel)) == 0){
  p_l_sel_p <- "no"
}else{
  p_l_sel_p <- ggarrange(plotlist = p_l_sel, ncol = 3, nrow = 2)
}



pdf(paste(args[3],args[4],"dabestr_high_sel.pdf",sep = ""),width = 14)
#pdf("Differential_abundant_ARGs/dabestr/6_overlap_args/ICU.dabestr_high_sel.pdf",width = 14)
p_h_sel_p
dev.off()



#pdf(paste(args[3],args[4],"dabestr_low.pdf",sep = ""),width = 14)
#p_l_p
#dev.off()


pdf(paste(args[3],args[4],"dabestr_low_sel.pdf",sep = ""),width = 14)
#message(names(p_l_sel_p))
#pdf("Differential_abundant_ARGs/dabestr/6_overlap_args/ICU.dabestr_low_sel.pdf",width = 14)
p_l_sel_p
dev.off()

#print(unpaired)

#paste(round(unpaired$result$difference,2),"\n", "[",unpaired$result$ci,"CI  ",round(unpaired$result$bca_ci_low,3),"; ",round(unpaired$result$bca_ci_high,3),"]",sep = "")



#p1 <- p_l[["STAPHYLOCOCCUS_MUPA_CONFERRING_RESISTANCE_TO_MUPIROCIN"]] 
#p2 <- p_l[["MULTIDRUG_ABC_TRANSPORTER"]]
#p3 <- p_l[["VANRI"]]
#p4 <- p_l[["VANZ"]] 
#p5 <- p_l[["TET40"]]
#p6 <- p_l[["MTRA"]]

#p <- ggarrange(p1,p2,p3,p4,p5,p6,widths = 10,ncol = 2,nrow=3)
#pdf("/Volumes/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/tidy_directory/script/Differential_abundant_ARGs/dabestr/6_overlap_args/6overlap.pdf")
#print(p)
#dev.off()

#save.image("/Volumes/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/tidy_directory/script/Differential_abundant_ARGs/dabestr/6_overlap_args/dabest.RData")



# ICU- vs Control 

a <- 0 
b <- 0

dabestr_list <- list()
dabestr_list_result <- list()
dabestr_list_result_sel <- list()
arg_sel <- vector()

p_h <- list() 
p_h_sel <- list() 
p_l <- list()
p_l_sel <- list()

j <- 1 
for (i in 1:length(ARGs_name)){
  arg <- ARGs_name[i]
  
  ARGs_RPKM_long_sel <- filter(ARGs_RPKM_long, ID == arg) 
  test_na <- filter(ARGs_RPKM_long_sel, Group %in% c("ICU-","Healthy"))
  
  if (mean(test_na$RPKM) == 0) {
    #ARGs_name <- ARGs_name[ARGs_name!= arg]
    next 
  }
  
  unpaired <- ARGs_RPKM_long_sel %>% 
    dabest(Group, RPKM, 
           idx = list(c("Healthy", "ICU-")),
           paired = FALSE
    )
  
  
  
  if (unpaired$result$difference > 0){
    a<- a+1 
    p_h[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"",sep = ""))
  }else if (unpaired$result$difference < 0){
    b <- b+1 
    p_l[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"",sep = "")) 
  }
  
  if ((unpaired$result$difference > 0 & unpaired$result$bca_ci_low >0) ||  (unpaired$result$difference < 0 & unpaired$result$bca_ci_high < 0)){
    arg_sel[j] <- arg
    dabestr_list_result_sel[[j]] <- unpaired$result[,c(1:13)]
    
    j <- j + 1 
  }
  if ((unpaired$result$difference > 0 & unpaired$result$bca_ci_low >0)){
    p_h_sel[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"",sep = "_"))
  }else if ((unpaired$result$difference < 0 & unpaired$result$bca_ci_high < 0)){
    p_l_sel[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"",sep = "_"))
  }
  
  dabestr_list[[i]] <- unpaired
  dabestr_list_result[[i]] <- unpaired$result
  
  
}



dabestr_out_sel <- do.call(rbind, dabestr_list_result_sel) 
dabestr_out_sel$arg <- arg_sel

a_sel <- do.call(rbind, dabestr_list_result_sel) 
a_sel$arg <- arg_sel

dim(dabestr_out_sel)

dabestr_out <- dabestr_out[order(dabestr_out$difference, decreasing = T),]


write_delim(dabestr_out_sel, paste(args[3],args[4],"dabeset.ICUmvsHealthy.output.txt",sep = ""), delim = "\t")



print(a)
print(b)

p_h_sel_p <- ggarrange(plotlist = p_h_sel, ncol = 3, nrow = 2)
p_l_sel_p <- ggarrange(plotlist = p_l_sel, ncol = 3, nrow = 2)


pdf(paste(args[3],args[4],"ICUm.dabestr_high_sel.pdf",sep = ""),width = 14)
#pdf("Differential_abundant_ARGs/dabestr/6_overlap_args/ICUmvsHealthy.dabestr_high_sel.pdf",width = 14)
p_h_sel_p
dev.off()


pdf(paste(args[3],args[4],"ICUm.dabestr_low_sel.pdf",sep = ""),width = 14)
#message(names(p_l_sel_p))
#pdf("Differential_abundant_ARGs/dabestr/6_overlap_args/ICUmvsHealthy.dabestr_low_sel.pdf",width = 14)
p_l_sel_p
dev.off()


# ICU+ vs Healthy 

a <- 0 
b <- 0

dabestr_list <- list()
dabestr_list_result <- list()
dabestr_list_result_sel <- list()
arg_sel <- vector()

p_h <- list() 
p_h_sel <- list() 
p_l <- list()
p_l_sel <- list()

j <- 1 
for (i in 1:length(ARGs_name)){
  arg <- ARGs_name[i]
  ARGs_RPKM_long_sel <- filter(ARGs_RPKM_long, ID == arg) 
  test_na <- filter(ARGs_RPKM_long_sel, Group %in% c("Healthy","ICU+"))
  
  if (mean(test_na$RPKM) == 0) {
    #ARGs_name <- ARGs_name[ARGs_name!= arg]
    next 
  }
  
  unpaired <- ARGs_RPKM_long_sel %>% 
    dabest(Group, RPKM, 
           idx = list(c("Healthy", "ICU+")),
           paired = FALSE
    )
  if (unpaired$result$difference > 0){
    a<- a+1 
    p_h[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"RPKM",sep = "_"))
  }else if (unpaired$result$difference < 0){
    b <- b+1 
    p_l[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"RPKM",sep = "_")) 
  }
  
  if ((unpaired$result$difference > 0 & unpaired$result$bca_ci_low >0) ||  (unpaired$result$difference < 0 & unpaired$result$bca_ci_high < 0)){
    arg_sel[j] <- arg
    dabestr_list_result_sel[[j]] <- unpaired$result[,c(1:13)]
    
    j <- j + 1 
  }
  if ((unpaired$result$difference > 0 & unpaired$result$bca_ci_low >0)){
    p_h_sel[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"RPKM",sep = "_"))
  }else if ((unpaired$result$difference < 0 & unpaired$result$bca_ci_high < 0)){
    p_l_sel[[arg]] <- plot(unpaired,group.summaries = TRUE, rawplot.ylabel = paste(arg,"RPKM",sep = "_"))
  }
  
  dabestr_list[[i]] <- unpaired
  dabestr_list_result[[i]] <- unpaired$result
  
  
}


dabestr_out_sel <- do.call(rbind, dabestr_list_result_sel) 
dabestr_out_sel$arg <- arg_sel

a_sel <- do.call(rbind, dabestr_list_result_sel) 
a_sel$arg <- arg_sel

dim(dabestr_out_sel)

dabestr_out <- dabestr_out[order(dabestr_out$difference, decreasing = T),]
#dim(dabestr_out_sel)

write_delim(dabestr_out_sel, paste(args[3],args[4],"dabeset.ICUpvsHealthy.output.txt",sep = ""), delim = "\t")



print(a)
print(b)

p_h_sel_p <- ggarrange(plotlist = p_h_sel, ncol = 3, nrow = 2)
p_l_sel_p <- ggarrange(plotlist = p_l_sel, ncol = 3, nrow = 2)


pdf(paste(args[3],args[4],"ICUp_dabestr_high_sel.pdf",sep = ""),width = 14)
#pdf("Differential_abundant_ARGs/dabestr/6_overlap_args/ICUpvsHealthy.dabestr_high_sel.pdf",width = 14)
p_h_sel_p
dev.off()


pdf(paste(args[3],args[4],"ICUp_dabestr_low_sel.pdf",sep = ""),width = 14)
#message(names(p_l_sel_p))
#pdf("Differential_abundant_ARGs/dabestr/6_overlap_args/ICUpvsHealthy.dabestr_low_sel.pdf",width = 14)
p_l_sel_p
dev.off()


save.image(paste(args[3],args[4],"dabestr.RData",sep = ""))
#save.image("dabestr.RData")


#load("/Volumes/Lucy/documents/201910_ICU_ITS/05_ARGsprofiling/deepARG/webserver_identity80_cov50/script/result/dabestr/202007/dabestr.Rdata")
