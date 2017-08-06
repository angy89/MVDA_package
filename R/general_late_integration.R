
#' A MVDA Function
#'
#' This function implements the general linear integration process
#' @param A is the matrix resulting from the concatenation of the membership matrix for each view
#' @param k the desidered number of cluster 
#' @param eps is the precision parameter
#' @param patients_classes is the factor of real class labels of the patients
#' @return a list containing the general linear integration results results, the clustering assignment, the confusion matrix, the error in confusion matrix, the matrix of clusters centroids, the label for each cluster and the matrix with fisher test results for each cluster
#' @export
#' 


general_late_integration <- function(A,k,alfa,eps,patientsDB,patients_classes){ 
  GLI(A=A,k=k,alfa=alfa,eps=eps)->GLI.res
  B <- GLI.res$B
  rownames(B) <- rownames(A)
  P <- GLI.res$P
  cluster <- unlist(apply(B,1,which.max))
  cm <- confusion_matrix(classes=patients_classes,clustering=cluster,patientsDB=t(patientsDB),nCluster=k)
  class_label(cm) -> ClassLabels
  
  error <- CMsup_error(cm,dim(A)[1])
  center <- findCenter(DB=t(patientsDB),clust_vector=cluster)
  
  fisher_test(cluster,ClassLabels,cm,patients_classes) -> fisher_test_res
  
  toRet <- list(GLI.res = GLI.res, cluster = cluster, cm = cm, error=error, center = center,labels = ClassLabels, fisher = fisher_test_res)
  return(toRet)
}