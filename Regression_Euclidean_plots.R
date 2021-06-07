#source('functions_needed.R')
#load("workspace for novels.RData")
source('PCA_Euclidean.R')

m<-dim(List_L[[1]])[1]
pca <- PCA_euc$pca

#times to make predictions at
times_dickens <- seq(1836,1870,1)
times_austen <- seq(1794, 1815,1)
ntms_D <- length(times_dickens)
ntms_A <- length(times_austen)

#____________________________________________________________________
######DICKENS#######
coeff_intercept_dic <- matrix(0,m,m)
coeff_diff_dic <- matrix(0,m,m)
for (i in 1:m){ 
  for (j in i:m){
    coeff_dic <- coefficients(lm(L[i,j,1:16] ~ times[1:16]))
    coeff_intercept_dic[i,j] <- coeff_dic[1]
    coeff_diff_dic[i,j] <- coeff_dic[2]
    coeff_intercept_dic[j,i] <- coeff_intercept_dic[i,j]
    coeff_diff_dic[j,i] <- coeff_diff_dic[i,j]
  }
}

#predicted values at times in 'times_dickens' are
predicted_glB_dic <- array(0, dim=c(m,m,ntms_D))
vectorized_pred_glB_dic <- matrix(0,m *m ,ntms_D)
for(i in 1:ntms_D){
  predicted_glB_dic[,,i] <- as.matrix(proj_sparse( coeff_intercept_dic+times_dickens[i]*coeff_diff_dic))
  vectorized_pred_glB_dic[,i] <- as.vector(predicted_glB_dic[,,i])
  }

#____________________________________________________________________
########AUSTEN############
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

#predicted values at times in 'times_austen' are
predicted_glB_aus <- array(0, dim=c(m,m,ntms_A))
vectorized_pred_glB_aus <- matrix(0,m *m ,ntms_A)
for(i in 1:ntms_A){
  predicted_glB_aus[,,i] <- as.matrix(proj_sparse( (coeff_intercept_aus+times_austen[i]*coeff_diff_aus)))
  vectorized_pred_glB_aus[,i] <- as.vector(predicted_glB_aus[,,i])
}

#____________________________________________________________________
###plots####
plot(pca$x[,1:2], col=1, ylab='coordinate 2', xlab='coordinate 1', cex=0.1)
text(pca$x[,1:2],labels=label, col=gray.colors(25, start = 0, end = 0.7)[rank(times)], cex=0.5)
abline(lm(predict(pca, newdata = t(vectorized_pred_glB_dic))[,2] ~ predict(pca, newdata = t(vectorized_pred_glB_dic))[,1])$coef, col=1, lty=2)
abline(lm(predict(pca, newdata = t(vectorized_pred_glB_aus))[,2] ~ predict(pca, newdata = t(vectorized_pred_glB_aus))[,1])$coef, col=1)
