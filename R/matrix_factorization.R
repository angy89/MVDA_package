#' A MVDA Function
#'
#' This function implements the matrix factorization integration process
#' @param A is the matrix resulting from the concatenation of the membership matrix for each view
#' @param k the desidered number of cluster
#' @param iter_mat the maximum allowed iteration
#' @param eps is the precision parameter
#' @param nView is the number of views
#' @param patientsDB is your matrix dataset where patients are represented by the prototypes selected in the previous step
#' @param patients_classes is the factor of real class labels of the patients
#' @return a list containing the matrix factorization results, the clustering assignment, the confusion matrix, the error in confusion matrix, the matrix T of views influence, the matrix of clusters centroids, the labels for each cluster and the matrix of fisher test for each cluster
#' @export
#'


MF <- function(k,A,eps,iter_max,nView,V1.lc1,V1.lc2,
                                 V1.lc3,V1.lc4,V1.lc5,patientsDB, patients_classes){
  matrix_factorization(X=t(A),k1=k,eps=eps,iter_max=iter_max)->mf
  cluster <- apply(mf$H,2,which.max)

  cm <- confusion_matrix(classes=patients_classes,clustering=cluster,patientsDB=t(patientsDB),nCluster=k)
  class_label(cm) -> ClassLabels
  error <- CMsup_error(cm,dim(A)[1])

  find_T(P=mf$P,k=k, nView=nView,V1.lc1=V1.lc1,V2.lc2=V1.lc2,V3.lc3=V1.lc3,V4.lc4 = V1.lc4,V5.lc5=V1.lc5) -> T_mat
  center <- findCenter(DB=t(patientsDB),clust_vector=cluster)
  fisher_test(cluster,ClassLabels,cm,patients_classes) -> fisher_test_res

  toRet <- list(mf_res = mf,mv_cluster=cluster,confMat = cm, error = error, T_mat = T_mat, center=center,labels = ClassLabels, fisher=fisher_test_res)
  return(toRet)
}
