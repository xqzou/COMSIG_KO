source("../MMRDetect/Gen_catalogues.R")
source("../MMRDetect/MMRDetect.classify.R")
source("../MMRDetect/MMRDetect.train.R")
source("../MMRDetect/MMRDetect.compute.variables.R")



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
