#source('functions_needed.R')
#load("workspace for novels.RData")

m<-dim(List_L[[1]])[1]

LC<-LC[1:m,1:m,]
words <- words[1:m]

#####################################################################
# test set up
#Dickens
DL <- L[,,1]*0
DL2 <- DL
CD <- LC[,,1]*0
for (i in 1:16){
  DL <- DL+((abs(L[,,i])))
  DL2 <- DL2+ ((abs(L[,,i]))**2)
  CD <- CD+abs(LC[,,i])
}
DL <- DL/16

DLvar <- (DL2 - DL^2*16)/16

#plot(DL,DLvar^0.5, ylab='sample standard deviation', xlab='sample mean', pch=19, cex=0.5)
#abline(c(0,0),0.2, col=2)
#linear regression downweight large values of mean relative to SD with five fold increase in mean)
srd <- DL*0.2
DLvar <- (184*srd**2 + 16*DLvar)/200 # 184*srd**2 is (1-w)(beta xbar)^2 part of equation



#Austen
AL <- L[,,1]*0
AL2 <- AL
CA <- LC[,,1]*0
for (i in 17:23){
  AL <- AL+((abs(L[,,i])))
  AL2 <- AL2+((L[,,i]**2))
  CA <- CA+abs(LC[,,i])
}
AL <- AL/7
ALvar <- (AL2 - AL^2*7)/7

#linear regression (by eye!)
psd <-  AL*0.2

#plot(AL,ALvar^0.5, ylab='sample standard deviation', xlab='sample mean', pch=19, cex=0.5)
#abline(c(0,0),0.2, col=2)
ALvar <- (193*psd**2+7*ALvar)/200
Poolsd<- sqrt( (16*DLvar + 7*ALvar)/21 )



##################################################
#make Poolsd huge for word combinations that do not appear in Austen at all 
#or in Dickens at all ,D since the variability is undefined

Poolsd[CD <= 0] <- 100
Poolsd[CA <= 0] <- 100


xi <- median(Poolsd[Poolsd<10])*sqrt(1/16+1/7)
print(xi)


#######################################
### Euclidean tests
T <- (abs(DL)-abs(AL))/(xi+Poolsd*sqrt(1/16+1/7) ) 

Teuclidean <- T
#matrix of test statistics

#no of pairs to display NP
NP <- 100


#####################################
#Significantly more in Dickens 
#get the top 100 pairs for drawing the graph
tem <- T-diag(diag(T))
Dthresh <- sort(tem[Poolsd<10],decreasing=TRUE)[NP*2]
print(Dthresh)

NT <- T*0
NT[T >= Dthresh ] <- 1
sum(NT)

A <-  NT - diag(diag(NT)) 

ones <- rep(1, dim(A)[1])
nzn <- (A%*%ones!=0)
AADickens <- A[nzn,nzn]

ggn <- network(AADickens)
set.seed(2)
plot(ggnet2(ggn,edge.alpha=0.5,label=words[nzn],color="yellow"))


#####################################
#significantly more in Austen
#get the top 100 pairs

tem <- T-diag(diag(T))
Athresh <- sort(tem[Poolsd<10],decreasing=FALSE)[NP*2]


print(Athresh)

NT <- NT*0
NT[T <= Athresh ] <- 1
sum(NT)

A <-  NT - diag(diag(NT)) 

nzn <- (A%*%ones!=0)
AAAusten <- A[nzn,nzn]

ggn <- network(AAAusten)
set.seed(2)
plot(ggnet2(ggn,edge.alpha=0.5,label=words[nzn],color="cyan"))
