#source('functions_needed.R')
#load("workspace for novels.RData")

m<-dim(List_L[[1]])[1]
Hel <- defh(m-1)


#____________________________________________________________________
#DICKENS

#_____________________________
#under null (coeff_diff_dic=D1==0)
D_0 <- matrix(0,m,m)
for (i in 1:m){
  for (j in 1:m){
    D_0[i,j] <- coefficients(lm(L[i,j,1:16] ~ 1))[1]
  }
}

diff_with_mean <- array(0, dim=c(m,m,16))
vec_diff_mean <- matrix(0,(m-1)*(m)/2, 16)
for (i in 1:16){
  diff_with_mean[,,i] <- L[,,i]-D_0
  vec_diff_mean[,i] <- vech(Hel%*%as.matrix(diff_with_mean[,,i])%*%t(Hel))
}

sigma_0 <- rep(0, (m-1)*(m)/2)
for (i in 1:((m-1)*(m)/2))
  sigma_0[i] <- (var(as.vector(vec_diff_mean[i,])))
cov_0 <- .sparseDiagonal(m*(m-1)/2,sigma_0    )
f_0 <- rep(0,16)
for (i in 1:16){
  f_0[i] <- (-1/2*t(vec_diff_mean[,i])%*%solve(cov_0)%*%vec_diff_mean[,i])-0.5*(sum(log(diag(cov_0))))-0.5*((m-1)*m/2)*log(2*pi)#as cov_0 diag eigenvalues are just diag
}

#_____________________________
#under alternative
coeff_intercept_dic <- matrix(0,m,m)
coeff_diff_dic <- matrix(0,m,m)
for (i in 1:m){
  for (j in i:m){
    coeff_intercept_dic[i,j] <- coefficients(lm(L[i,j,1:16] ~ times[1:16]))[1]
    coeff_diff_dic[i,j] <- coefficients(lm(L[i,j,1:16] ~ times[1:16]))[2]
    coeff_intercept_dic[j,i] <- coeff_intercept_dic[i,j]
    coeff_diff_dic[j,i] <- coeff_diff_dic[i,j]
  }
}

diff_with_pred <- array(0, dim=c(m,m,16))
vec_diff_pred <- matrix(0, (m-1)*(m)/2, 16)
for (i in 1:16){
  diff_with_pred[,,i] <- L[,,i]-coeff_intercept_dic-times[i]*coeff_diff_dic
  vec_diff_pred[,i] <- vech(Hel%*%as.matrix(diff_with_pred[,,i])%*%t(Hel))
  }

sigma_1 <- rep(0, (m-1)*(m)/2)
for (i in 1:((m-1)*(m)/2))
  sigma_1[i] <- (var(as.vector(vec_diff_pred[i,])))
cov_1 <- .sparseDiagonal(m*(m-1)/2,sigma_1    )
f_1 <- rep(0,16)
for (i in 1:16)
  f_1[i] <- (-1/2*t(vec_diff_pred[,i])%*%solve(cov_1)%*%vec_diff_pred[,i])-0.5*(sum(log(diag(cov_1))))-0.5*((m-1)*m/2)*log(2*pi)#exclude (2*pi)^(m^2)* as will dive by this


#_____________________________
#TEST STATISTIC
-2*(sum(f_0)-sum(f_1))
pchisq(-2*(sum(f_0)-sum(f_1)), m*(m-1)/2) #1-this is p-val
p_val_dick_euc <- 1-pchisq(-2*(sum(f_0)-sum(f_1)), m*(m-1)/2) #if this is small then evidence that regression exists
print(c('the p value for existance of linear regression for the Dickens sample is ', p_val_dick_euc))
qchisq(0.95,m*(m-1)/2)




#____________________________________________________________________
#Austen

#_____________________________
#under null (coeff_diff_dic=D1==0)
D_0 <- matrix(0,m,m)
for (i in 1:m){
  for (j in 1:m){
    D_0[i,j] <- coefficients(lm(L[i,j,17:23] ~ 1))[1]
  }
}


diff_with_mean <- array(0, dim=c(m,m,7))
vec_diff_mean <- matrix(0,(m-1)*(m)/2, 7)
for (i in 1:7){
  diff_with_mean[,,i] <- L[,,i+16]-D_0
  vec_diff_mean[,i] <- vech(Hel%*%as.matrix(diff_with_mean[,,i])%*%t(Hel))
}

sigma_0 <- rep(0, (m-1)*(m)/2)
for (i in 1:((m-1)*(m)/2))
  sigma_0[i] <- (var(as.vector(vec_diff_mean[i,])))
cov_0<-.sparseDiagonal(m*(m-1)/2,sigma_0    )
f_0<-rep(0,7)
for (i in 1:7)
  f_0[i] <- (-1/2*t(vec_diff_mean[,i])%*%solve(cov_0)%*%vec_diff_mean[,i])-0.5*(sum(log(diag(cov_0))))#-0.5*((m-1)*m/2)*log(2*pi)can exclude as constant


#_____________________________
#under alternative
coeff_intercept_aus <- matrix(0,m,m)
coeff_diff_aus <- matrix(0,m,m)
for (i in 1:m){
  for (j in i:m){
    coeff_intercept_aus[i,j] <- coefficients(lm(L[i,j,17:23] ~ times[17:23]))[1]
    coeff_diff_aus[i,j] <- coefficients(lm(L[i,j,17:23] ~ times[17:23]))[2]
    coeff_intercept_aus[j,i] <- coeff_intercept_aus[i,j]
    coeff_diff_aus[j,i] <- coeff_diff_aus[i,j]
  }
}

diff_with_pred <- array(0, dim=c(m,m,7))
vec_diff_pred <- matrix(0, (m-1)*(m)/2, 7)
for (i in 1:7){
  diff_with_pred[,,i] <- L[,,i+16]-coeff_intercept_aus-times[i+16]*coeff_diff_aus
  vec_diff_pred[,i] <- vech(Hel%*%as.matrix(diff_with_pred[,,i])%*%t(Hel))
  }

sigma_1<-rep(0, (m-1)*(m)/2)
for (i in 1:((m-1)*(m)/2))
  sigma_1[i] <- (var(as.vector(vec_diff_pred[i,])))
cov_1<-.sparseDiagonal(m*(m-1)/2,sigma_1    )
f_1<-rep(0,7)
for (i in 1:7)
  f_1[i] <- (-1/2*t(vec_diff_pred[,i])%*%solve(cov_1)%*%vec_diff_pred[,i])-0.5*(sum(log(diag(cov_1))))#-0.5*((m-1)*m/2)*log(2*pi)can exclude as constant


#_____________________________
#TEST STATISTIC
-2*(sum(f_0)-sum(f_1))
pchisq(-2*(sum(f_0)-sum(f_1)), m*(m-1)/2)
p_val_aust_euc <- 1-pchisq(-2*(sum(f_0)-sum(f_1)), m*(m-1)/2)#if this is small then evidence that regression exists
print(c('the p value for existance of linear regression for the Austen sample is ', p_val_aust_euc))
qchisq(0.95,m*(m-1)/2)


