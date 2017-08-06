#' A MVDA Function
#'
#' This function evaluate the VAL index for different clustering algorithms and divverent value of K for the number of clusters
#' @param DB your matrix dataset, where the rows are the features and the columns are the patients
#' @param nCentersRange is a numeric vector of numbers of clusters
#' @param clustAlgos is a vector with the names of the clustering algorithms. Possible values are KMeans,Pam,Ward,Pvclust,SOM
#' @param som_dim is a list with couples of dimensions for the som e.g. list(c(2,5),c(2,6))
#' @return a dataframe with information on the clustering results
#' @export
best_preprocessing = function(DB,nCentersRange,clustAlgos,som_dim ){
  
  optAlgos = c("KMeans","Pam","Ward","Pvclust","SOM")
  id_check = clustAlgos %in% optAlgos
  
  if(sum(id_check) < length(clustAlgos)){
    print(paste(clustAlgos[which(id_check==F)], "is not in the list of available clustering algorithms!\n"))
    print("Choose one of the following algorithms: ")
    cat(optAlgos)
    return(-1)
  }
  
  if(length(nCentersRange) != length(som_dim)){
    print("The list of K values differs in length with the list of SOM dimensions")
    return(-1)
  }
  
  clSummaryDat = c()
  
  for(i in 1:length(nCentersRange)){
    
    for(clAlgo in clustAlgos){
      if(clAlgo == "KMeans"){
        prep_res = kmeans_preprocessing(DB = DB,nCenters = nCentersRange[i])
      }
      if(clAlgo == "Pam"){
        prep_res = pamk_preprocessing(DB = DB,nCenters = nCentersRange[i])
      }
      if(clAlgo ==" Ward"){
        prep_res = ward_hls_prep(DB = DB,nCenters = nCentersRange[i])
      }
      if(clAlgo == "Pvclust"){
        prep_res = pvclust_prep(DB = t(DB),nCenters = nCentersRange[i])
      }
      if(clAlgo == "SOM"){
        prep_res = SOM_preprocessing(DB = DB,xdim = som_dim[[i]][1],ydim = som_dim[[i]][2])
      }
      
      cl_summary <- clustering_summary(DB=DB, cluster=prep_res$clustering)
      clSummaryDat = rbind(clSummaryDat,c(clAlgo,cl_summary))
    }  
  }
  
  clSummaryDat = clSummaryDat[order(clSummaryDat[,2]),]
  colnames(clSummaryDat) = c("Algo","nClusters",colnames(cl_summary)[2:6])
  clSummaryDat =  as.data.frame(clSummaryDat)
  clSummaryDat$Index = round(as.numeric(as.vector(clSummaryDat$Index)),2)
  
  g = ggplot(clSummaryDat, aes(Algo, Index, group =nClusters,colour =nClusters)) + 
    geom_line() 

  return(list(clSummaryDat=clSummaryDat,gplt = g))
  
}