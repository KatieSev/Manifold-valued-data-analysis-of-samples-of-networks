#source('functions_needed.R')
#load("workspace for novels.RData")

m<-dim(List_L[[1]])[1]

PCA_euc <- PCA_GL(List_L, TRUE,FALSE, FALSE)
plot(PCA_euc$pca$x[,1:2], col=author_col, ylab='coordinate 2', xlab='coordinate 1', cex=0.1, xaxp  = c(-0.02, 0.02, 2), yaxp  = c(-0.01, 0.01, 2))
text(PCA_euc$pca$x[,1:2],labels=label, col=rainbow(25)[rank(times)], cex=0.5) #WAS 0.75 save 260

plot(c(cumsum(PCA_euc$pca$sdev^2)/sum(PCA_euc$pca$sdev^2)), ylim=c(0,1), xlab='PC number', ylab='Cumulative % variance explained', pch=19, cex=0.75)


#####

dickens_eucmean <- Mean_GL(List_L[1:16], euc=TRUE, sqrt=FALSE, proc=FALSE)
austen_eucmean <- Mean_GL (List_L[17:23], euc=TRUE, sqrt=FALSE, proc=FALSE)

##################
#Barplot to look at degree in PCs
PC1_euc <-matrix(PCA_euc$pca$rotation[,1],m,m)
rownames(PC1_euc)<-lookup[1:m]
colnames(PC1_euc)<-lookup[1:m]

order_weight_euc<-rbind(sort(diag(PC1_euc)/sum(diag(abs(PC1_euc)))),
                        diag(austen_eucmean-dickens_eucmean)[order(diag(PC1_euc)/sum(diag(abs(PC1_euc))))]/sum(diag(abs(austen_eucmean-dickens_eucmean))))

barpos<-barplot(as.vector(order_weight_euc)[c(1:20, (2*m-19):(2*m))],beside=T,col=c('black','gray'), names.arg = rep("",40),
                ylim=c(-0.1,0.09))
axis(1, lwd.tick=0, labels=FALSE)

for (i in 1:10){
  axis(1, at = barpos[2*i-1], label="", tck = -0.01)
  axis(1, at=barpos[2*i-1],  lwd = 0,labels=lookup[order(diag(PC1_euc)/sum(diag(abs(PC1_euc))))][i],las=3, cex.axis=0.6, line = 0.5)
}

for (i in (m-9):m){
  axis(1, at = barpos[2*i-2*(m-20)-1], label="", tck = -0.01)
  axis(1, at=barpos[2*i-2*(m-20)-1],lwd=0, line=0.5, labels=lookup[order(diag(PC1_euc)/sum(diag(abs(PC1_euc))))][i],las=2, cex.axis=0.6)
  
}


##################
#Barplot to look at degree in PCs
PC2_euc <-matrix(PCA_euc$pca$rotation[,2],m,m)
rownames(PC2_euc)<-lookup[1:m]
colnames(PC2_euc)<-lookup[1:m]

order_weight_euc_2<-rbind(sort(diag(PC2_euc)/sum(diag(abs(PC2_euc)))))

barpos<-barplot(as.vector(order_weight_euc_2)[c(1:10, (m-9):m)],beside=T,col=c('black'), names.arg = rep("",20),
                ylim=c(-0.15,0.09))
axis(1, lwd.tick=0, labels=FALSE)

for (i in 1:10){
  axis(1, at = barpos[i], label="", tck = -0.01)
  axis(1, at=barpos[i],  lwd = 0,labels=lookup[order(diag(PC2_euc)/sum(diag(abs(PC2_euc))))][i],las=3, cex.axis=0.6, line = 0.5)
}

for (i in (m-9):m){
  axis(1, at = barpos[i-(m-20)], label="", tck = -0.01)
  axis(1, at=barpos[i-(m-20)],lwd=0, line=0.5, labels=lookup[order(diag(PC2_euc)/sum(diag(abs(PC2_euc))))][i],las=2, cex.axis=0.6)
  
}
