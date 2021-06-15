source('functions_needed.R')
load("workspace for novels.RData")

topwords <- 100 #code in paper is for topwords=1000
L<- L[1:topwords,1:topwords, ]
for (i in 1:23){
  List_L[[i]] <- List_L[[i]][1:topwords,1:topwords]
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

