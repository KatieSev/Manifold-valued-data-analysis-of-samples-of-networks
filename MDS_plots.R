load("workspace for novels.RData")
library(shapes)
nt<-23

Dmat<-matrix(0,nt,nt)
for (i in 1:(nt-1)){
  print(i)
  for (j in (i+1):nt){
    Dmat[i,j]<-sqrt(sum((L[,,i]-L[,,j])^2))
    Dmat[j,i]<-Dmat[i,j]
  }
}
plot(hclust(as.dist(Dmat),"ward.D2"))
plot(cmdscale(Dmat),type="n")
text(cmdscale(Dmat),labels=1:nt)
lines(cmdscale(Dmat))
Dmateucl<-Dmat



nt<-nr
Dmat<-matrix(0,nt,nt)
for (i in 1:(nt-1)){
  print(i)
  for (j in (i+1):nt){
    Dmat[i,j]<-norm( rootmat(L[,,i])- rootmat(L[,,j]) , type='f')
    Dmat[j,i]<-Dmat[i,j]
  }
}
plot(hclust(as.dist(Dmat),"ward.D2"))
plot(cmdscale(Dmat),type="n")
text(cmdscale(Dmat),labels=1:nt)
lines(cmdscale(Dmat))
Dmatsqrt<-Dmat

nt<-nr
Dmat<-matrix(0,nt,nt)
for (i in 1:(nt-1)){
  print(i)
  for (j in (i+1):nt){
    Dmat[i,j]<-distcov( L[,,i],L[,,j] , method="Procrustes")
    Dmat[j,i]<-Dmat[i,j]
  }
}
plot(hclust(as.dist(Dmat),"ward.D2"))
plot(cmdscale(Dmat),type="n")
text(cmdscale(Dmat),labels=1:nt)
lines(cmdscale(Dmat))
Dmatproc<-Dmat


#save as 230 by 230 pixels
novcol<-c(rep(2,times=16),rep(4,times=7))
dickens_austen_ids[17:23]<-c("em","pe","pr","ls","ma","no","se")

dend<-hclust(as.dist(Dmateucl),method="ward.D2") 
plot(dend,labels=dickens_austen_ids, xlab="", main="", sub="", cex=0.5, xaxt='n')
abline(h=0.08, lty=2)
set.seed(1)
plot(cmdscale(Dmateucl),type="n",ylab="MDS2", xaxt='n', yaxt='n', xlab="")
text(cmdscale(Dmateucl),labels=dickens_austen_ids,col=novcol, cex=0.5, xaxt='n', yaxt='n')


dend<-hclust(as.dist(Dmatsqrt),method="ward.D2") 
plot(dend,labels=dickens_austen_ids, xlab="", main="", sub="", cex=0.5)
abline(h=0.7, lty=2)
set.seed(1)
plot(cmdscale(Dmatsqrt),type="n",xlab="",ylab="MDS2",yaxp  = c(-0.1, 0.1, 2), xaxt='n', yaxt='n')
text(cmdscale(Dmatsqrt),labels=dickens_austen_ids,col=novcol, cex=0.5, xaxt='n', yaxt='n')

set.seed(1)
dend<-hclust(as.dist(Dmatproc),method="ward.D2") 
plot(dend,labels=dickens_austen_ids, xlab="", main="", sub="", cex=0.5)
abline(h=0.7, lty=2)
plot(cmdscale(Dmatproc),type="n",xlab="MDS1",ylab="MDS2",yaxp  = c(-0.1, 0.1, 2), xaxt='n', yaxt='n')
text(cmdscale(Dmatproc),labels=dickens_austen_ids,col=novcol, cex=0.5, xaxt='n', yaxt='n')
