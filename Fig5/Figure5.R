source("../MMRDetect/Gen_catalogues.R")
source("../MMRDetect/MMRDetect.classify.R")
source("../MMRDetect/MMRDetect.train.R")
source("../MMRDetect/MMRDetect.compute.variables.R")
library(pROC)
library(ROSE)

col_ff_MMRDetect_var <- read.table("Colorectal_train_test.txt", sep = "\t", header = T, as.is = T)
rownames(col_ff_MMRDetect_var) <- col_ff_MMRDetect_var$AnonymousID
col_ff_MMRDetect_var_train <- col_ff_MMRDetect_var[col_ff_MMRDetect_var$dataset=="train", c("RepIndel_num", "Del_rep_mean", "MMR_sum", "maxcossim","MSI_status")] 
col_ff_MMRDetect_var_train$MSI_status <- factor(col_ff_MMRDetect_var_train$MSI_status)

col_ff_MMRDetect_var_test <- col_ff_MMRDetect_var[col_ff_MMRDetect_var$dataset=="test", c("RepIndel_num", "Del_rep_mean", "MMR_sum", "maxcossim","MSI_status")] 
col_ff_MMRDetect_var_test$MSI_status <- factor(col_ff_MMRDetect_var_test$MSI_status)

set.seed(123)
# Normalize
col_ff_MMRDetect_var_train2 <- col_ff_MMRDetect_var_train
col_ff_MMRDetect_var_train2$RepIndel_num <- col_ff_MMRDetect_var_train2$RepIndel_num/max(col_ff_MMRDetect_var_train2$RepIndel_num)
col_ff_MMRDetect_var_train2$MMR_sum <- col_ff_MMRDetect_var_train2$MMR_sum/max(col_ff_MMRDetect_var_train2$MMR_sum)
glm_model<- glm(MSI_status~.,data=col_ff_MMRDetect_var_train2,family = binomial(link = "logit"))
saveRDS(glm_model,"./MMRDetect.rds")
#col_ff_MMRDetect_var_test2 <- col_ff_MMRDetect_var_test
#col_ff_MMRDetect_var_test2$RepIndel_num <- col_ff_MMRDetect_var_test2$RepIndel_num/max(col_ff_MMRDetect_var_test2$RepIndel_num)
#col_ff_MMRDetect_var_test2$MMR_sum <- col_ff_MMRDetect_var_test2$MMR_sum/max(col_ff_MMRDetect_var_test2$MMR_sum)
col_ff_MMRDetect_var_test_classified <- MMRDetect.classify(col_ff_MMRDetect_var_test,glm_model)
write.table(col_ff_MMRDetect_var_test_classified, "glm_test_col.txt", sep = "\t", col.names = T, row.names = T, quote = F)

col_ff_MMRDetect_var_train_classified <- MMRDetect.classify(col_ff_MMRDetect_var_train,glm_model)
write.table(col_ff_MMRDetect_var_train_classified, "glm_train_col.txt", sep = "\t", col.names = T, row.names = T, quote = F)

col_ff_MMRDetect_var_classified <- rbind(col_ff_MMRDetect_var_test_classified, col_ff_MMRDetect_var_train_classified)
col_ff_MMRDetect_var_classified$Sample <- rownames(col_ff_MMRDetect_var_classified)
write.table(col_ff_MMRDetect_var_classified, "glm_train_test_col.txt", sep = "\t", col.names = T, row.names = T, quote = F)


###############################################################
# Test with corrected MSI status ROC of IHC, MMRDetect, MSIseq
###############################################################
col_3results <- read.table("./Results/Colorectal_3results_anonymous.txt", sep = "\t", header = T, as.is = T)
col_3results$real_status <- col_3results$MMRDetect
pdf(file="roc_curve_3methods.pdf", onefile=TRUE,height=5,width=5, useDingbats=FALSE)
ROSE::roc.curve(col_3results$real_status, col_3results$MSI_status)
ROSE::roc.curve(col_3results$real_status, col_3results$MMRDetect,add.roc=T, col=2, lwd=2, lty=2)
ROSE::roc.curve(col_3results$real_status, col_3results$MSIseq, add.roc=T, col=3, lwd=3, lty=3)
legend("bottomright", c("IHC (AUC:0.979)", "MMRDetect (AUC:1)", "MSIseq (AUC:0.996)"), col=1:3, lty=1:3, lwd=2)
dev.off()

###############################################################
# Test with corrected MSI status ROC of IHC, MMRDetect, MSIseq
###############################################################



