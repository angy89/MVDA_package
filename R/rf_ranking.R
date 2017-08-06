#' A MVDA Function
#'
#' This function determines a ranking of predictors by using random foresto from the package randomForest.
#' @param prototype is the matrix of prototype, patients are on columns and featuers are on rows
#' @param patient_classes is the factor of class label
#' @param MTRY vector of possible values for the random forest's mtry parameters
#' @param NTREE vector of possible values for the random forest's ntree parameters
#' @param NPERM vector of possible values for the random forest's nperm parameters
#' @param NODESIZE vector of possible values for the random forest's nodesize parameters
#' @keywords multi-view clustering; feature ranking; randomForest;
#' @return A matrix of importance measure, one row for each predictor variable. The column(s) are different importance measures.
#' @export
#'
#'
rf_ranking <- function(prototype,patient_classes, MTRY=c(1,2,3),NTREE=c(100,200,500),
                       NPERM = c(1,2,3), NODESIZE=c(1,5)){
  require(randomForest)
  require(rminer)

  if(is.null(prototype) || is.null(patient_classes)){
    stop("You must insert Prototypes and Class labels!\n")
  }else{
  N_patient <- dim(prototype)[2]

patients <- t(prototype)
N_FEAT = dim(patients)[2]

patients = cbind(patients,patient_classes)

colnames(patients) = c(colnames(patients)[1:N_FEAT],"breast")
gsub("_","UNDESC",colnames(patients)) -> colnames(patients)
gsub("-","_",colnames(patients)) -> colnames(patients)
gsub("[*]","_STAR",colnames(patients)) -> colnames(patients)
gsub("[.]","_DOT",colnames(patients)) -> colnames(patients)
gsub("[:]","_COLON",colnames(patients)) -> colnames(patients)

patients = as.data.frame(patients)
patients$breast = as.factor(patients$breast)

search=list(smethod="grid",search=list(mtry=MTRY,ntree=NTREE,nPerm = NPERM, nodesize=NODESIZE),
            convex=0,metric="AUC",method=c("kfold",3,12345))
M=fit(breast~., as.data.frame(patients),model="randomForest",search=search,fdebug=TRUE)
print(M@mpar)
M@object$importance->imp
rf2 = M@object

gsub("_STAR","*",rownames(imp)) -> rownames(imp)
gsub("_DOT",".",rownames(imp)) -> rownames(imp)
gsub("_COLON",":",rownames(imp)) -> rownames(imp)
gsub("_","-",rownames(imp)) -> rownames(imp)
gsub("UNDESC","_",rownames(imp)) -> rownames(imp)
return(imp)
  }
}
