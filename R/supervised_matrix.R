#' A MVDA Function
#'
#' This function create the membership matrix of patient and classes
#' @param patient_classes is the vector with class labels
#' @param nRow is the number of rows of the matrix (number of patients)
#' @param nCol is the number of columns of the matrix (number of classes)
#' @param row_id are the rows names
#' @param col_id are the cols names
#' @keywords multi-view clustering; confusion matrix
#' @return a membership matrix
#' @export
#' 

supervised_matrix <- function(classes,nRow,nCol,row_id,col_id){
  
  SUP <- matrix(0,nrow=nRow,ncol=nCol) 
  for(i in 1:nRow){
    SUP[i,classes[i]]<-1
  }
  rownames(SUP)<- row_id
  colnames(SUP)<- col_id
  return(SUP)
}