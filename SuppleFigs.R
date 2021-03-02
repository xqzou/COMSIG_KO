# Supplementary Figures

# Perform 10 fold cross validation
col_ff_MMRDetect_var_train3 <- col_ff_MMRDetect_var_train2[sample(nrow(col_ff_MMRDetect_var_train2)),]
# Create 10 equally size folds
folds <- cut(seq(1,nrow(col_ff_MMRDetect_var_train3)),breaks = 10, labels = F)
coef_10 <- NULL
for(i in 1:10){
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- col_ff_MMRDetect_var_train3[testIndexes, ]
  trainData <- col_ff_MMRDetect_var_train3[-testIndexes, ]
  set.seed(123)
  glm_model_i <- glm(MSI_status~., data=trainData, family = binomial)
  glm_prob_i <- predict.glm(glm_model_i, testData[,-which(names(testData)=="MSI_status")], type="response")
  test_roc <- roc(testData$MSI_status ~glm_prob_i, plot=F, print.auc=F)
  print(paste0(i,"_AUC:", as.numeric(test_roc$auc)))
  print(glm_model_i$coefficients)
  coef_10 <- rbind(coef_10, c(glm_model_i$coefficients, as.numeric(test_roc$auc)))
}
coef_10 <- as.data.frame(coef_10)
names(coef_10) <- c("intercept","RepIndel_num","Del_rep_mean","MMR_sum","maxcossim","auc")
write.table(coef_10,"coef_10.txt",sep = "\t", col.names = T, row.names = F, quote = F)
coef_10_melt <- melt(coef_10[,-c(1,6)])

pdf(file="variable_cv10_v2.pdf", onefile=TRUE,width = 6,height = 3, useDingbats = F)
g <-ggplot(coef_10_melt, aes(y=variable, x=value, fill=variable)) + geom_boxplot()
g <- g+xlab("Coefficients")+ylab("Variable")+xlim(-60,0)
g <- g+theme(axis.text.x=element_text(size=10,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
print(g)
#ggMarginal(g, type = "histogram", fill="transparent")
dev.off()

