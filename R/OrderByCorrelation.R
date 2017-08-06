#' A MVDA Function
#'
#' This function allows you to rank feature relevance for each multi-view cluster by correlation.
#' It works by removing a feature at a time and calculating correlation between patients in the cluster. 
#' The more the correlation decrease, the more the feature is relevant for the clusters.
#' It also order patients by correlation.
#' @param DB your matrix dataset
#' @param cluster.vector a numeric vector of clustering results
#' @keywords multi-view clustering; feature-relevance; correlation
#' @return a list with two object sorted_features is the list of feature ordered by correlation and sorted_patients that is the list of patients ordered by correlations
#' @export

OrderByCorrelationFeaturesPatients <- function(DB, cluster.vector){
    orderedPatientByCorrelation <- list()
    for(i in names(table(cluster.vector))){
        index <- which(cluster.vector==i)
        if(length(index)==1){
            patient_corr <- 1;
            names(patient_corr) <- colnames(DB)[index]
            orderedPatientByCorrelation <- lappend(orderedPatientByCorrelation,patient_corr)
        }else{
            elemClust <- DB[,index]
            matCor <- cor(elemClust, method="pearson");
            patient_corr <- (apply(matCor,1,FUN=function(riga){
                sum(riga)
            }))
            patient_corr <-sort(patient_corr)
            orderedPatientByCorrelation <- lappend(orderedPatientByCorrelation,patient_corr)
        }
    }
    
    orderedFeaturesByCorrelation <- list();
    clustCorr <- c()
    #se quando tolgo una feature la correlazione diminuisce allora la feature era altamente significativa
    for(i in names(table(cluster.vector))){ #per ogni cluster
        #cat("cluster1,",i,"\n")
        index <- which(cluster.vector==i);
        if(length(index)==1){
            clustCorr <- c(clustCorr,1)
        }else{
            clustCorr <- c(clustCorr,mean(cor(DB[,index],method="pearson")))
        }
        featCor <-NULL
        for(j in 1:dim(DB)[1]){ #per ogni feature
            #cat("feature ",j,"\n")
            feat <- 1:dim(DB)[1];
            feat <- feat[-j]
            which(cluster.vector==i)->index;
            if(length(index)==1){
                featCor <- c(featCor,1)
            }else{
                elemClust <- DB[feat,index]
                matCor <- cor(elemClust, method="pearson");
                featCor <- c(featCor,mean(matCor))
            }
        }
        names(featCor) <- rownames(DB)    
        featCor = sort(featCor,decreasing = TRUE)
        orderedFeaturesByCorrelation <- lappend(orderedFeaturesByCorrelation,featCor)
    }
    
    toRet <- list(sorted_features= orderedFeaturesByCorrelation,sorted_patients = orderedPatientByCorrelation,corr_all_prototypes=clustCorr)
    return(toRet)
}
