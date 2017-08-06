#' A MVDA Function
#'
#' This function assign labels to each cluster. the label is assigned to the majority. That is, for each cluster identifies the predominance of objects of a given class. The cluster will be assigned the same label of the predominant class.
#' 
#' @param CM is the confusion matrix
#' @return a vector containing a label for each cluster
#' @export


class_label <- function(CM){
  CM <- as.matrix(CM)
  cluster_class_labels <- apply(CM,1,which.max)
  return(cluster_class_labels)
}


