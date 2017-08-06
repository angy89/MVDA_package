#' A MVDA Function
#'
#' This function evaluate the error in a confusion matrix
#' @param CM_sup confusion matrix
#' @param nPat number of patients
#' @keywords multi-view clustering; confusion matrix
#' @return an error value
#' @export
CMsup_error <- function(CM_sup,nPat){
  apply(CM_sup,1,FUN=function(row){
    sum(row[-which.max(row)])
  })->somma
  
  return(sum(somma)/nPat)
}