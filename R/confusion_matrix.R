#' A MVDA Function
#'
#' This function construct confusion matrix between patient classes and the obtained clustering
#' @param classes is a vector of patient labels
#' @param clustering is a vector of clustering results
#' @param patientsDB is your matrix dataset
#' @param nCluster is the number of obtained clusters
#' @keywords multi-view clustering; confusion matrix
#' @return the confusion matrix 
#' @export

confusion_matrix <- function(classes, clustering, patientsDB, nCluster){
  
  gene_clust <- clustering
  nClass <- length(table(classes))
  
  matrix_summary <-matrix(classes,ncol=1,nrow = length(classes));
  rownames(matrix_summary)<- rownames(patientsDB);
  
  summary <- matrix_summary
  for(i in 1:length(table(classes))){
    index <- which(matrix_summary[,1]==attr(table(classes)[i],which="name"));
    matrix_summary[index,1] <-i;
  }
  
  ConfusionMatrix_gene <- matrix(0,nrow=nCluster,ncol=nClass)
  for(i in 1:nCluster){
    tab.i <- table(matrix_summary[which(gene_clust==i),1]);
    for(j in 1:length(tab.i)){
      ConfusionMatrix_gene[i,as.integer(attr(tab.i[j],which="name"))]<-tab.i[j]
    }
  }
  rownames(ConfusionMatrix_gene) <- paste("cluster",1:nCluster,sep="")
  colnames(ConfusionMatrix_gene) <- paste("classi",1:nClass,sep="")
  
  return(ConfusionMatrix_gene)
}