findPvclustCenter <- function(DB,cluster.pv,n_cluster){
	center <- NULL #matrice dei centroidi
	name_center <- NULL
	for(i in 1:n_cluster){
		clustElem <- DB[,cluster.pv$clusters[[i]]];
        matCor <- cor(clustElem, method="pearson");
        bestCenter <- which.max(apply(matCor,1,FUN=function(riga){
            sum(riga)
        }))
		center <- rbind(center,clustElem[,bestCenter])
		name_center <- c(name_center,attr(bestCenter,which="names"))
	}
	rownames(center)<-name_center
	return(center)
}


findVector <- function(DB,cluster.pv){
	cluster.vector <- rep(0,length(colnames(DB)))
	for(i in  1:length(colnames(DB))){
  		for(j in 1:length(cluster.pv$clusters)){
    		if(colnames(DB)[i] %in% cluster.pv$clusters[[j]]){
      		cluster.vector[i] <- j
    		}
  		}
  	}
	return(cluster.vector)
}

pvcluster <- function(DB,hclust.method,nboot,dist.method,alpha,show.print,r,printLocation){
  cat("calculate Pvcluster\n")
	genes.pv <- pvclust(DB, method.hclust=hclust.method, nboot=nboot, method.dist=dist.method,r=r)
	if(show.print){
		plot(genes.pv)
	}
  cat("cat dendrogram\n")
	cluster.pv <- pvpick(genes.pv, alpha=alpha)
	n_cluster <- length(cluster.pv$clusters)

	cluster.vector <- findVector(DB,cluster.pv)  	
  	center <- findPvclustCenter(DB,cluster.pv,n_cluster)
    
    #SALVO IL RISULTATO DEL PVCLUST SEMPLICE
    pvclust.cluster.pv <- cluster.pv
    
    #creo la lista che memorizzi quali geni sono del vecchio risultato del pvclust e quasi sono quelli che aggiungo nel passo EM
    typeGene <- list();
    for(i in 1:n_cluster){
        zero_v <- rep(1,length(pvclust.cluster.pv$clusters[[i]]))
        typeGene <- lappend(typeGene,zero_v)
    }
    
	#EFFETTUO UN PASSO EM SU I GENI NON CLUSTERIZZATI PER TROVARE IL CENTROIDE PIù VICINO 
	#(con correlazione più alta) AL QUALE ASSEGNARLO ED INOLTRE AGGIORNO I CENTRI DEI NUOVI CLUSTER
	#cenco di assegnare i geni non clusterizzati ai cluster più vicini
	#la vicinanza è calcolata tramite la correlazione
  
  if(length(which(cluster.vector==0))>0){
      cat("EM step\n")
    	notClustered <-DB[,which(cluster.vector==0)]
    	resumeBestCenter <- NULL
    	for(i in 1:dim(notClustered)[2]){
    		bestCenter <- which.max(apply(center,1,FUN=function(elem){
    			cor(elem,notClustered[,i])
    		}))	
    		resumeBestCenter <- c(resumeBestCenter,bestCenter)
    		cluster.pv$clusters[[bestCenter]]<-c(cluster.pv$clusters[[bestCenter]],colnames(notClustered)[i])
            typeGene[[bestCenter]]<-c(typeGene[[bestCenter]],3)
    	}
        old_center <- center
    	center <- findPvclustCenter(DB,cluster.pv,n_cluster)
    	cluster.vector <- findVector(DB,cluster.pv)
      toRet <- list(pvclustRes=genes.pv,n_cluster = n_cluster,typeGene = typeGene,old_center=old_center,pvclust.cluster.pv=pvclust.cluster.pv,cluster.pv = cluster.pv, cluster.vector = cluster.vector, center = center,genes.pv=genes.pv)
      
  }
  else{
    toRet <- list(pvclustRes=genes.pv,n_cluster = n_cluster,typeGene = typeGene,old_center=center,pvclust.cluster.pv=pvclust.cluster.pv,cluster.pv = cluster.pv, cluster.vector = cluster.vector, center = center,genes.pv=genes.pv)
  }
	tab <- table(cluster.vector)
    cat("End pvcluster \n")
	return(toRet)
}

#ORDINO I GENI E LE FEATURE DI OGNI CLUSTER IN BASE ALLA CORRELAZIONE
OrderByCorrelation <- function(n_cluster, cluster.pv, DB, cluster.vector){
	orderedGenesByCorrelation <- list()
	OrderedVector <- list();
	for(i in 1:n_cluster){
		elemClust <- DB[,which(cluster.vector==i)]
		matCor <- cor(elemClust, method="pearson");
        gene_corr <- (apply(matCor,1,FUN=function(riga){
            sum(riga)
        }))
        gene_corr <-sort(gene_corr)
        orderedGenesByCorrelation <- lappend(orderedGenesByCorrelation,gene_corr)
        
	}
	
	orderedFeaturesByCorrelation <- list();
	#se quando tolgo una feature la correlazione diminuisce allora la feature era altamente significativa
	for(i in 1:n_cluster){ #per ogni cluster
		
		featCor <-NULL
		for(j in 1:dim(DB)[1]){ #per ogni feature
			feat <- 1:dim(DB)[1];
			feat <- feat[-j]
			elemClust <- DB[feat,which(cluster.vector==i)]
			matCor <- cor(elemClust, method="pearson");
        	featCor <- c(featCor,mean(matCor))
		}
		
		names(featCor) <- rownames(DB)
		featCor <- sort(featCor)
		
		orderedFeaturesByCorrelation <- lappend(orderedFeaturesByCorrelation,featCor)
	}
	
	toRet <- list(feature= orderedFeaturesByCorrelation,genes = orderedGenesByCorrelation)
}


printClusterDifference <- function(DB,cluster.pv,center,pvclust.cluster.pv,old_center,directory_name,n_cluster,typeGene,centroidLwd){
    for(i in 1:n_cluster){
        
        par(mfrow=c(1,2))
        file_name <- paste(directory_name,i,".pdf",sep="");
        pdf(file = file_name);
        #FIGURA 1
        myCol = typeGene[[i]][1]
        plot(DB[,pvclust.cluster.pv$clusters[[i]][1]],type="l", col=myCol,ylim=c(3,15), main=paste("Previous pvclust result: Cluster: ",i))
        
        for(j in 2:length(pvclust.cluster.pv$clusters[[i]])){
            myCol = typeGene[[i]][j]
            lines(DB[,pvclust.cluster.pv$clusters[[i]][j]], col=myCol)
        }
        
        #calcolo il centroide come la media dei valori
        lines(old_center[i,],col="red", lwd= centroidLwd,type="o")
        #legend(x = 10,y = 14,legend=paste("Numbers of elemnts in the cluster: ",length(cluster.pv$clusters[[i]])))
        legend(x = "topright",legend=c("Previous Results","EM-step","Centroid"),fill=c(1,3,"red"))
        
        #FIGURA 2
        myCol = typeGene[[i]][1]
        plot(DB[,cluster.pv$clusters[[i]][1]],type="l", col=myCol,ylim=c(3,15), main=paste("Modified pvclust result: Cluster: ",i))
        
        for(j in 2:length(cluster.pv$clusters[[i]])){
            myCol = typeGene[[i]][j]
            lines(DB[,cluster.pv$clusters[[i]][j]], col=myCol)
        }
        
        #calcolo il centroide come la media dei valori
        lines(center[i,],col="red", lwd=centroidLwd,type="o")
        #legend(x = 10,y = 14,legend=paste("Numbers of elemnts in the cluster: ",length(cluster.pv$clusters[[i]])))
        legend(x = "topright",legend=c("Previous Results","EM-step","Centroid"),fill=c(1,3,"red"))
        dev.off();
    }
    
}

#DB è considerato come non trasposto
adjacency_matrix <- function(DB,cluster.vector,n_cluster){
	n_elem <- length(cluster.vector)
	adjmat <- matrix(data = 0,nrow = n_elem,ncol = n_cluster);
	rownames(adjmat)<- rownames(DB); #il nome degli oggetti
	for(i in 1:n_elem){
		adjmat[i,cluster.vector[i]]<-1
	}	
	
	return(adjmat)
}




