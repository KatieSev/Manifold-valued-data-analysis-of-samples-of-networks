#source('functions_needed.R')
#load("workspace for novels.RData")

m<-dim(List_L[[1]])[1]

nt<-23

novcol<-c(rep(2,times=16),rep(4,times=7))


#Euclidean
Dmat<-matrix(0,nt,nt)
for (i in 1:(nt-1)){
  print(i)
  for (j in (i+1):nt){
    Dmat[i,j]<-sqrt(sum((L[,,i]-L[,,j])^2))
    Dmat[j,i]<-Dmat[i,j]
  }
}
Dmateucl<-Dmat
dend<-hclust(as.dist(Dmateucl),method="ward.D2") 
plot(dend,labels=label, xlab="", main="", sub="", cex=0.5, xaxt='n')
abline(h=0.08, lty=2)
set.seed(1)
plot(cmdscale(Dmateucl),type="n",ylab="MDS2", xaxt='n', yaxt='n', xlab="")
text(cmdscale(Dmateucl),labels=label,col=novcol, cex=0.5, xaxt='n', yaxt='n')




#Squareroot Euclidean
Dmat<-matrix(0,nt,nt)
L_sqrt <- L
for (i in 1:23){
  L_sqrt[,,i] <- rootmat(L[,,i])
}
for (i in 1:(nt-1)){
  print(i)
  for (j in (i+1):nt){
    Dmat[i,j]<-norm( L_sqrt[,,i] - L_sqrt[,,j] , type='f')
    Dmat[j,i]<-Dmat[i,j]
  }
}
Dmatsqrt<-Dmat
dend<-hclust(as.dist(Dmatsqrt),method="ward.D2") 
plot(dend,labels=label, xlab="", main="", sub="", cex=0.5)
abline(h=0.7, lty=2)
set.seed(1)
plot(cmdscale(Dmatsqrt),type="n",xlab="",ylab="MDS2",yaxp  = c(-0.1, 0.1, 2), xaxt='n', yaxt='n')
text(cmdscale(Dmatsqrt),labels=label,col=novcol, cex=0.5, xaxt='n', yaxt='n')




#Procrustes size-and-shape
Dmat<-matrix(0,nt,nt)
for (i in 1:(nt-1)){
  print(i)
  for (j in (i+1):nt){
    Dmat[i,j]<-distcov( L[,,i],L[,,j] , method="Procrustes")
    Dmat[j,i]<-Dmat[i,j]
  }
}
Dmatproc<-Dmat
dend<-hclust(as.dist(Dmatproc),method="ward.D2") 
plot(dend,labels=label, xlab="", main="", sub="", cex=0.5)
abline(h=0.7, lty=2)
set.seed(1)
plot(cmdscale(Dmatproc),type="n",xlab="MDS1",ylab="MDS2",yaxp  = c(-0.1, 0.1, 2), xaxt='n', yaxt='n')
text(cmdscale(Dmatproc),labels=label,col=novcol, cex=0.5, xaxt='n', yaxt='n')
