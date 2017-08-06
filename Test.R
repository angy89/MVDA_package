#setwd("~/Dropbox/paper tesi/bmc_article/Script_paper/MVDA/")

library(MVDA)
data(tcga_breast)

# 1 step miRNASeq Preprocessing
miRNA_km <- kmeans_preprocessing(DB = miRNA, nCenters = 10)

prep_res = best_preprocessing(DB = miRNA,nCentersRange = c(5,10,20,30,40),
                              som_dim = list(c(1,5),c(2,5),c(4,5),c(6,5),c(5,8)),
                              clustAlgos = c("KMeans","Pam","Ward","SOM"))
plot(prep_res$gplt)

miRNA_km_summary <- clustering_summary(DB=miRNA, cluster=miRNA_km$clustering)


rf_rank_res = rf_ranking(prototype=miRNA_km$centers,patient_classes=info$x, 
                         MTRY=c(1,2,3),NTREE=c(100,200,500), NPERM = c(1,2,3), NODESIZE=c(1,5))

RNA_pamk <- pamk_preprocessing(DB = RNA,nCenters = 50)
RNA_pamk_summary <- clustering_summary(DB=RNA, cluster=RNA_pamk$clustering)

# 2 step prototype ranking
sda_km_miRNA <- sda_ranking(prototype = miRNA_km$centers,
                            classes = info$x,ranking.score = "avg")
plot(sda_km_miRNA)
sda_RNA <- sda_ranking(prototype = RNA_pamk$centers,
                       classes = info$x,ranking.score = "avg")
plot(sda_RNA)

#3 step single view clustering
prototypes_miRNA <- miRNA_km$centers[rownames(sda_km_miRNA),]
nCenter <- 4
miRNA_pat_km <- kmeans_sv(nCenters = nCenter,prototype = prototypes_miRNA)
cm_km <- confusion_matrix(classes = info$x,
	clustering = miRNA_pat_km$clustering,
	patientsDB = t(prototypes_miRNA),
	nCluster = nCenter)
miRNAerror <- CMsup_error(CM_sup = cm_km,nPat = dim(prototypes_miRNA)[2])
prototypes_RNA <- RNA_pamk$centers[rownames(sda_RNA),]
nCenter <- 4
RNA_pat_km <- kmeans_sv(nCenters = nCenter,prototype = prototypes_RNA)
cm_km <- confusion_matrix(classes = info$x,
						  clustering = RNA_pat_km$clustering,
						  patientsDB = t(prototypes_RNA),
						  nCluster = nCenter)
RNAerror <- CMsup_error(CM_sup = cm_km,nPat = dim(prototypes_RNA)[2])

#4 step integration
nrows <- dim(info)[1]
ncols <- length(table(miRNA_pat_km$clustering))
rows_id <- rownames(info)
cols_id <- paste("Clustering",names(table(miRNA_pat_km$clustering)),sep=" ")
X1 <- supervised_matrix(classes = miRNA_pat_km$clustering,nRow = nrows,
						nCol = ncols,row_id = rows_id,col_id  = cols_id)
X2 <- supervised_matrix(classes = RNA_pat_km$clustering,nRow = nrows,
						nCol = ncols,row_id = rows_id,col_id  = cols_id)
SUP <- supervised_matrix(classes = as.numeric(info$x),nRow = nrows,
						nCol = length(table(info$x)),
						row_id = rows_id,col_id  = cols_id)
X <- cbind(X1,X2,SUP)
K <- dim(X)[2]

patientsDB <- t(cbind(t(prototypes_RNA),t(prototypes_miRNA)))
mf_res <- MF(A = X,k = 10,eps = 0.000001,iter_max = 1000,nView = 3,V1.lc1 = 4,
			V1.lc2 = 4,V1.lc3 = 4,V1.lc4 = 0,V1.lc5 = 0,
			patients_classes = info$x,patientsDB = patientsDB)
GLI_res <- general_late_integration(A = X,k = K,alfa = 1,
 			eps = 0.001,patientsDB = patientsDB,patients_classes = info$x)

plotCM(mf_res$confMat)

feature_importance_for_clusters = OrderByCorrelationFeaturesPatients(DB=patientsDB,cluster.vector = mf_res$mv_cluster)#,n_cluster=nrow(mf_res$confMat))
stability_results_MF = mvclst(method = "MF",nPatients = 151,nClusters = K,A = X,eps = 0.000001,iter_max = 1000)
#stability_results_GLI = mvclst(method = "GLI",A = X,nPatients = 151,nClusters = K,eps = 0.001,alfa =1)


var_col = plot_feature_variance_between_centroids(center=mf_res$center,v_limit = 12,h_limit=1)


center_correlation(centers = mf_res$center)

mv_clust_boxplot(centers = mf_res$center, 
                 clusterMV = mf_res$mv_cluster,
                 cluster_class= mf_res$labels,
                 var_col=var_col,nFeat = 6,patientsDB=patientsDB,
                              patClasses=c("Basal","Her2","LumA","LumB","Normal"),
                              patClassColors = c("purple","blue","green","red","orange"))
