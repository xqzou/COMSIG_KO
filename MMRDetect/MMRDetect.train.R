# Train MMRDetect



MMRDetect.train <- function(mutationVariable, classification, cancerType = NULL) {
  
  
  ## match the data with classification
  trainset = mutationVariable
  trainset = merge(trainset, classification[,c("Sample","MSI_status")], by="Sample")
  
  if(nrow(trainset)<50){
    warning('training set size < 50')  
  }
  
  # normalize RepIndel_num and MMR_sum
  trainset$RepIndel_num <- trainset$RepIndel_num/max(trainset$RepIndel_num)
  trainset$MMR_sum <- trainset$RepIndel_num/max(trainset$MMR_sum)
  
  ## build model with trainset
  trainset$MSI_status<-as.factor(trainset$MSI_status)
  glm_model_logit = glm(MSI_status~., data = trainset, family = binomial(link="logit"))
  glm_model_logit
  
}

