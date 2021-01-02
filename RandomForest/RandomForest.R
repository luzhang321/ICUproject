# Random Forest model 
# classfier icu- vs icu+ using top20 features & check how many present in previous sig.pathway analysis 
# name : cuiqin 

library(pacman)
p_load(tidyverse,
       data.table,
       caret,
       skimr,
       plotROC,
       mlbench,
       MLeval,
       ggpubr,
       pROC,
       PRROC)

# working directory 
# ----------------------------------------------------------------------------------------
setwd("/media/lu/Lucy/documents/201910_ICU_ITS/10_machinelearning/data/")

# read file 
# ----------------------------------------------------------------------------------------
all_data_combine <- fread("bacteria_pathway_arg_fungi_group.txt")

# filter data 
#----------------------------------------------------------------------------------------------------------------------------

all_data_combine_filter_no_intersection <- select(all_data_combine,  !matches("\\|g__")) %>%
  select(., !matches("\\|unclassified")) # 1995 

# only select bacteria & pathway information 


all_data_combine_filter_no_intersection <- all_data_combine_filter_no_intersection[,c(1:913,1995)]

# replace ICU- with ICUm and ICU+ with ICUp 
all_data_combine_filter_no_intersection[all_data_combine_filter_no_intersection$Bgroup == "ICU+",]$Bgroup <- "ICUp"
all_data_combine_filter_no_intersection[all_data_combine_filter_no_intersection$Bgroup == "ICU-",]$Bgroup <- "ICUm"


all_data_input_ICUpICUm_no_intersection <- filter(all_data_combine_filter_no_intersection, Bgroup == "ICUp" | Bgroup == "ICUm") %>%
  column_to_rownames(.,var="sample") 


# set train parameter 
ctrl <- trainControl(method="cv", 
                     number = 5,
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)

# train model : ICUp and ICUm 
set.seed(1234)
rf_ICUp_ICUm_model_no_intersection <- 
  train(Bgroup ~ ., 
        data=all_data_input_ICUpICUm_no_intersection, 
        method="rf", preProc=c("center", "scale","nzv"), 
        trControl=ctrl,
        tuneLength = 10,
        metric = "ROC"
  )


feature_importance_no_intersection <- varImp(rf_ICUp_ICUm_model_no_intersection)


importance_feature_no_intersection <- data.frame(feature_importance_no_intersection$importance,feature = rownames(feature_importance_no_intersection$importance))
# select top20 features 
importance_feature_top20_no_intersection <- head(importance_feature_no_intersection[order(importance_feature_no_intersection[,1],decreasing = T),],n=20) %>%
  .$feature %>%
  gsub("`", '', .)  

# train model again using this top20 features 
rf_ICUp_ICUm_data_top20_no_intersection <- all_data_input_ICUpICUm_no_intersection[,importance_feature_top20_no_intersection] 
rf_ICUp_ICUm_data_top20_no_intersection$Bgroup <- all_data_input_ICUpICUm_no_intersection$Bgroup

set.seed(1234)
rf_ICUp_ICUm_model_top20_no_intersection <- train(Bgroup ~ ., 
                                  data=rf_ICUp_ICUm_data_top20_no_intersection, 
                                  method="rf", preProc=c("center", "scale","nzv"), 
                                  trControl=ctrl,
                                  tuneLength = 10,
                                  metric = "ROC"
)



rf_ICUp_ICUm_model_top20_no_intersection$finalModel # this one you get an accuracy 
a_no_intersection <- filter(rf_ICUp_ICUm_model_top20_no_intersection$pred, mtry ==4)
a_no_intersection_confusion <- confusionMatrix(a_no_intersection$pred,a_no_intersection$obs,positive = "ICUp")
evalm(rf_ICUp_ICUm_model_top20_no_intersection)

# draw confusion matrix 
draw_confusion_matrix <- function(cm) {
  
  
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n', bty ="none")
  #title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  # rect(150, 430, 240, 370, col='#3F97D0')
  rect(150, 430, 240, 370, col='#F7AD50')
  text(195, 435, 'ICU-', cex=1.8)
  rect(250, 430, 340, 370, col='white')
  text(295, 435, 'ICU+', cex=1.8)
  text(125, 370, 'Prediction', cex=1.3, srt=90, font=2)
  text(245, 450, 'Reference', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='white')
  #rect(250, 305, 340, 365, col='#3F97D0')
  rect(250, 305, 340, 365, col='#F7AD50')
  text(140, 400, 'ICU-', cex=1.8, srt=90)
  text(140, 335, 'ICU+', cex=1.8, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='black')
  text(195, 335, res[2], cex=1.6, font=2, col='black')
  text(295, 400, res[3], cex=1.6, font=2, col='black')
  text(295, 335, res[4], cex=1.6, font=2, col='black')
  
  
}  




# rf_ICUp_ICUm_model_no_intersection$finalModel
# b <- filter(rf_ICUp_ICUm_model_no_intersection$pred, mtry == 911)
# confusionMatrix(b$pred,b$obs,positive = "ICUp")

pdf("../output/15.final.ROC.pdf")
rf_ICUp_ICUm_model_ROC <- roc(a_no_intersection$obs, a_no_intersection$ICUp) %>%
  plot(.,print.auc=TRUE,print.thres=FALSE, legacy.axes = FALSE)
dev.off()


library(PRROC)
library(ROCR)
data(ROCR.simple)
df <- data.frame(ROCR.simple)
pred <- prediction(df$predictions, df$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)

# PRROC_obj <- roc.curve(scores.class0 = df$predictions, weights.class0=df$labels,
#                        curve=TRUE)

a_no_intersection$test <- ifelse(a_no_intersection$obs == "ICUp",1,0 )
PRROC_obj <- PRROC::roc.curve(scores.class0 = a_no_intersection$ICUp, weights.class0=a_no_intersection$test,
                       curve=TRUE,rand.compute = TRUE)

pdf("../output/15.final.confusionmatrix.pdf")
draw_confusion_matrix(a_no_intersection_confusion)
dev.off()
pdf("../output/15.final.color.ROC.pdf",width = 8)
par(cex.axis=2, cex.lab =1.6,font.axis =1.6)
plot(PRROC_obj,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
text(1, 0.5, paste('AUC =', ' ',round(PRROC_obj$auc, digits = 3), sep =""), cex=1.4, font=2, col='black')
dev.off()
#predict(rf_ICUp_ICUm_model_top20_no_intersection,test_data_Bgroup_with_path)
varImp(rf_ICUp_ICUm_model_top20_no_intersection)





p <- plot(varImp(rf_ICUp_ICUm_model_top20_no_intersection, scale = TRUE))

rf_ICUp_ICUm_model_top20_no_intersection$finalModel
pdf("../output/15.final.top20.features.pdf",width = 10)
p
dev.off()


# features that present in the significant pathways from previous analysis 
# -----------------------------------------------------------------------------
sig_pathways <- readxl::read_xlsx("significant_list/significant_pathways_wilcoxon_full.xlsx")
sig_path_ICUm_ICUp <- sig_pathways$`ICU+ vs ICU-`[1:46]

common_sig_path <- intersect(importance_feature_top20_no_intersection,sig_path_ICUm_ICUp) 

setdiff(importance_feature_top20_no_intersection,sig_path_ICUm_ICUp)
# pathways classfication 

pathway_classification <- readxl::read_xlsx("significant_list/significant_pathways_wilcoxon_full_classification.xlsx")


common_sig_path[common_sig_path %in% pathway_classification$aminoacids] #10 
common_sig_path[common_sig_path %in% pathway_classification$SCFA] # 3 
common_sig_path[common_sig_path %in% pathway_classification$BA] #2 

#maunally paste it to a excel 
#

# how many in AAs, SCFAs and BAs pathways 

importance_feature_top20_no_intersection[importance_feature_top20_no_intersection %in% pathway_classification$aminoacids] #10
importance_feature_top20_no_intersection[importance_feature_top20_no_intersection %in% pathway_classification$SCFA] # 4
importance_feature_top20_no_intersection[importance_feature_top20_no_intersection %in% pathway_classification$BA] #2 


# significant species from previous analysis 

sig_species <- readxl::read_xlsx("significant_list/wilcox_diff_species.xlsx")
sig_species$`ICU- - ICU+`


save.image("/media/lu/Lucy/documents/201910_ICU_ITS/10_machinelearning/15_independent_corhort.RData")
load("/media/lu/Lucy/documents/201910_ICU_ITS/10_machinelearning/15_independent_corhort.RData")
