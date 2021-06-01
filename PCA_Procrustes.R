source('functions_needed.R')
load("workspace for novels.RData")


#PERFORM PCA
pca_proc <- PCA_GL(List_L, euc=FALSE, sqrt=FALSE, proc=TRUE)

plot(pca_proc$pca$x[,1], pca_proc$pca$x[,2], col=author_col, ylab='coordinate 2', xlab='coordinate 1', cex=0.1, xaxp  = c(-0.08, 0.16, 3), yaxp  = c(-0.08, 0.08, 2))
text(pca_proc$pca$x[,1], pca_proc$pca$x[,2],labels=label, col=rainbow(25)[rank(times)], cex=0.5) 

plot(c(0,cumsum(pca_proc$pca$sdev^2)/sum(pca_proc$pca$sdev^2)), ylim=c(0,1), xlab='PC number', ylab='Cumulative % variance explained', pch=19, cex=0.75)
