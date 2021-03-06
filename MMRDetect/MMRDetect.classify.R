#' MMRDetect classifier
#'
#' @param mutationVariable A list of input variables,"Del_rep_mean","RepIndel_num","MMR_sum","maxcossim"
#' @param classifier provided classifier
#' @return classification result

#' @export
MMRDetect.classify <- function(mutationVariable, classifier = MMRDclassifier) {
  classifyset = mutationVariable[,c("Del_rep_mean","RepIndel_num","MMR_sum","maxcossim")]
  
  classifyset$RepIndel_num <- classifyset$RepIndel_num/max(classifyset$RepIndel_num)
  classifyset$MMR_sum <- classifyset$MMR_sum/max(classifyset$MMR_sum)
  
  mutationVariable$glm_prob = predict.glm(classifier, newdata=classifyset, type="response")
  return(mutationVariable)
} 
