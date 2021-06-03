source('PCA_Squareroot.R')

m <- 1000

#times to make predictions at
times_dickens <- seq(1836,1870,1)
times_austen <- seq(1794, 1815,1)
ntms_D <- length(times_dickens)
ntms_A <- length(times_austen)


#set up
gl_sqrt <- array(0, dim=c(m,m,23))
for(i in 1:23)
  gl_sqrt[,,i] <- rootmat(L[,,i])

sqrtmean_all <- (Mean_GL(List_L, euc=FALSE, sqrt=TRUE, proc=FALSE))
sqrted_sqrtmean_all <- rootmat(sqrtmean_all)

#____________________________________________________________________
######DICKENS#######
coeff_intercept_dic <- matrix(0,m,m)
coeff_diff_dic <- matrix(0,m,m)

for (i in 1:m) {
  for (j in i:m) {
    coeff_intercept_dic[i, j] <- coefficients(lm(gl_sqrt[i, j, 1:16] ~ times[1:16]))[1]
    coeff_diff_dic[i, j] <- coefficients(lm(gl_sqrt[i, j, 1:16] ~ times[1:16]))[2]
    coeff_intercept_dic[j, i] <- coeff_intercept_dic[i, j]
    coeff_diff_dic[j, i] <- coeff_diff_dic[i, j]
  }
}
#predicted values at times in 'times_dickens' are
predicted_gl_sq_B_dic <- array(0, dim = c(m, m, ntms_D))
vectorized_pred_gl_sqB_dic <- matrix(0, m * m , ntms_D)
for (i in 1:ntms_D){
  predicted_gl_sq_B_dic[, , i] <- as.matrix(proj_sparse(( G_2(coeff_intercept_dic + times_dickens[i] * coeff_diff_dic)
                                                          %*%t(G_2(coeff_intercept_dic + times_dickens[i] * coeff_diff_dic)) )))
  vectorized_pred_gl_sqB_dic[, i] <-  as.vector(rootmat(predicted_gl_sq_B_dic[, , i])-sqrted_sqrtmean_all)
}

#____________________________________________________________________
######AUSTEN#######
coeff_intercept_aus <- matrix(0, m, m)
coeff_diff_aus <- matrix(0, m, m)
for (i in 1:m){
  for (j in i:m){
    coeff_intercept_aus[i, j] <- coefficients(lm(gl_sqrt[i, j, 17:23] ~ times[17:23]))[1]
    coeff_diff_aus[i, j] <- coefficients(lm(gl_sqrt[i, j, 17:23] ~ times[17:23]))[2]
    coeff_intercept_aus[j, i] <- coeff_intercept_aus[i, j]
    coeff_diff_aus[j, i] <- coeff_diff_aus[i, j]
  }
}

#predicted values at times in 'times_austen' are
predicted_gl_sq_B_aus <- array(0, dim=c(m,m,ntms_A))
vectorized_pred_gl_sq_B_aus <- matrix(0,m *m ,ntms_A)
for(i in 1:ntms_A){
  predicted_gl_sq_B_aus[,,i] <- as.matrix(proj_sparse((G_2(coeff_intercept_aus+times_austen[i]*coeff_diff_aus)
                                                       %*%t(G_2(coeff_intercept_aus+times_austen[i]*coeff_diff_aus))) ))
  vectorized_pred_gl_sq_B_aus[,i] <- as.vector(rootmat(predicted_gl_sq_B_aus[,,i])-sqrted_sqrtmean_all)
}

#____________________________________________________________________
######PLOT#######
plot(pca_sqrt$pca$x[,1:2], col=1, ylab='coordinate 2', xlab='coordinate 1', cex=0.1)
text(pca_sqrt$pca$x[,1:2],labels=label, col=gray.colors(25, start = 0, end = 0.7)[rank(times)], cex=0.5)
abline(lm(predict(pca_sqrt$pca, newdata = t(vectorized_pred_gl_sqB_dic))[,2] ~ predict(pca_sqrt$pca, newdata = t(vectorized_pred_gl_sqB_dic))[,1])$coef, col=1, lty=2)
abline(lm(predict(pca_sqrt$pca, newdata = t(vectorized_pred_gl_sq_B_aus))[,2] ~ predict(pca_sqrt$pca, newdata = t(vectorized_pred_gl_sq_B_aus))[,1])$coef, col=1)
