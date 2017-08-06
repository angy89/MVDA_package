#' A MVDA Function
#'
#' This function determines a ranking of predictors by computing CAT scores (correlation-adjusted t-scores) between the group centroids and the pooled mean by using fucntion of the package sda.
#' @param prototype is the matrix of prototype
#' @param classes is the factor of class label
#' @param fdr compute FDR values and HC scores for each feature. default = FALSE
#' @param ranking.score how to compute the summary score for each variable from the CAT scores of all classes. default = "avg"
#' @param lambda Shrinkage intensity for the correlation matrix. If not specified it is estimated from the data. lambda=0 implies no shrinkage and lambda=1 complete shrinkage. default = 0.5
#' @keywords multi-view clustering; feature ranking; cat-score;
#' @return return a matrix with the following columns: idx    original feature number  score  sum of the squared CAT scores across groups - this determines the overall ranking of a feature cat for each group and feature the cat score of the centroid versus the pooled mean
#' @export
#' 
sda_ranking <- function(prototype=NULL,classes=NULL,fdr=FALSE,ranking.score="avg",lambda=0.5){
  library(sda)

  if(is.null(prototype) || is.null(classes)){
    stop("You must insert Prototypes and Class labels!\n")
  }else{
    sda.ranking(Xtrain=t(prototype),fdr=FALSE,ranking.score="avg",L=classes,lambda=0.5) -> prototype.sda.rank
    return(prototype.sda.rank)
  }
}

