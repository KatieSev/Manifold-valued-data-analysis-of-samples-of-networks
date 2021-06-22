source('functions_needed.R')
load("workspace for novels.RData")


topwords <- 100 #code in paper is for topwords=1000
LC <- array(0, dim=c(topwords, topwords, 23))
L <- array(0, dim=c(topwords, topwords, 23))
#below converts LC and L to truncated version based on topwords value chosen
for (i in 1:23){
  List_LC[[i]] <- as.matrix(List_LC_full[[i]][1:topwords,1:topwords])
  List_L[[i]] <- as.matrix(List_LC_full[[i]][1:topwords,1:topwords])
  diag(List_L[[i]]) <- 0
  diag(List_L[[i]]) <- -rowSums(List_L[[i]])
  List_L[[i]] <- List_L[[i]]/sum(diag(List_L[[i]]))
  #
  LC[,,i] <- List_LC[[i]]
  L[,,i] <- List_L[[i]]
}



#MDS and dendrogram plots
source('MDS_plots.R')

#creating means, can then be saved and plots produced in cytoscape
source('Means.R')

#interpolation plots
source('Interpolation_plots.R')

#PCA plots
source('PCA_Euclidean.R')
source('PCA_Squareroot.R')
source('PCA_Procrustes.R')

#PCA with regression line added
source('Regression_Euclidean_plots.R')
source('Regression_Squareroot_plots.R')
source('Regression_Procrustes_plots.R')

#Regression Wilks test
source('Regression_Wilks_Theorem.R')

#Comparing means with z-type test
source('comparing_Austen_Dickens_means.R')

#Permutation test
source('permutation_test_Austen_Dickens.R')

