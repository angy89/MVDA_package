#' A MVDA Function
#'
#' run the fisher's exact test for each cluster in order to assign statistical significance to the label that has been assigned.
#' @param CM_sup is the confusion matrix
#' @param clust is the vector of clustering assignment
#' @param cl.cl is the vector of labels for each cluster
#' @param patients_classes is the factor of patient classes
#' @return a matrix containing a pvalue for each cluster for each class label.
#' @export
#' 
fisher_test <- function(clust,cl.cl,CM_sup, patients_classes){ 
  require(MASS)
  
  nCluster <- dim(CM_sup)[1]
  nClass <- length(names(table(patients_classes)))
  
  colSums(CM_sup) -> elem_par_class
  fisherTest <- matrix(0,nrow=nCluster,ncol=nClass)
  
  for(i in 1:nCluster){
    for(j in 1:nClass){
      cl1.test <- matrix(c(CM_sup[i,j],sum(CM_sup[i,-j]),elem_par_class[j],sum(elem_par_class[-j])),nrow=2,ncol=2,byrow=TRUE)
      fisherTest[i,j]<- fisher.test(cl1.test)$p.value
    }  
  }
  
  cbind(fisherTest,cl.cl)
  
  cl.cl <- names(table(patients_classes))[cl.cl]
  fisherTest <- cbind(fisherTest,cl.cl)
  colnames(fisherTest) <- c(names(table(patients_classes)),"Cluster class")
  rownames(fisherTest) <- paste("Cluster",1:nCluster,sep="")
  
  return(fisherTest)
}

