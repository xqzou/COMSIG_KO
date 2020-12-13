MMRDetect.classify <- function(mutationVariable, classifier = MMRDclassifier) {
  classifyset = mutationVariable[,c("Del_rep_mean","RepIndel_num","MMR_sum","maxcosim")]
  
  
  classifyset$glm_prob = predict.glm(classifier, newdata=classifyset, type="response")
  return(classifyset)
} 
