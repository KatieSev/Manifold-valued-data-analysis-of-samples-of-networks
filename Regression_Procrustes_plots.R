source('PCA_Procrustes.R')

#Result 2.2.2 in thesis is important on why we dont need helmet matrix
m <- 1000
times_dickens <- seq(1836,1870,1)
times_austen <- seq(1794, 1815,1)
ntms_D <- length(times_dickens)
ntms_A <- length(times_austen)

#set up
gl_proc_dic <- array(0, dim=c(m,m,23))
gl_proc_aus <- array(0, dim=c(m,m,23))

gl_proc_dic_1 <- estSS3(L[,,1:16])$tan
gl_proc_aus_1 <- estSS3(L[,,17:23])$tan

for (i in 1:16){
  gl_proc_dic[,,i] <- as.matrix(gl_proc_dic_1[,i], m,m)
}
for (i in 1:7){
  gl_proc_aus[,,i] <- as.matrix(gl_proc_aus_1[,i], m,m)
}

gl_proc_dic_pole <- estSS3(L[,,1:16])$mshape
gl_proc_aus_pole <- estSS3(L[,,17:23])$mshape

procmean_all <- (Mean_GL(List_L, euc=FALSE, sqrt=FALSE, proc=TRUE))
sqrted_procmean_all <- rootmat(procmean_all)

#____________________________________________________________________
######DICKENS#######
coeff_intercept_dic <- matrix(0,m,m)
coeff_diff_dic <- matrix(0,m,m)

for (i in 1:m) {
  for (j in i:m) {
    coeff_intercept_dic[i, j] <- coefficients(lm(gl_proc_dic[i, j, 1:16] ~ times[1:16]))[1]
    coeff_diff_dic[i, j] <- coefficients(lm(gl_proc_dic[i, j, 1:16] ~ times[1:16]))[2]
    coeff_intercept_dic[j, i] <- coeff_intercept_dic[i, j]
    coeff_diff_dic[j, i] <- coeff_diff_dic[i, j]
  }
}

#predicted values at times in 'times_dickens' are
predicted_gl_sq_B_dic <- array(0, dim = c(m, m, ntms_D))
vectorized_pred_gl_sqB_dic <- matrix(0, m * m , ntms_D)
for (i in 1:ntms_D){
  dummy <- coeff_intercept_dic + times_dickens[i] * coeff_diff_dic+rootmat(gl_proc_dic_pole)
  matched_dummy <- procOPA(rootmat(gl_proc_dic_pole), dummy, scale = FALSE, reflect = TRUE)$Bhat
  predicted_gl_sq_B_dic[, , i] <- as.matrix( proj_sparse( t(matched_dummy) %*% (matched_dummy) ) )
  vectorized_pred_gl_sqB_dic[, i] <-  as.vector(procOPA(sqrted_procmean_all,
                                                        rootmat(predicted_gl_sq_B_dic[, , i]),
                                                        scale=FALSE, reflect=TRUE)$Bhat - sqrted_procmean_all)
}


#____________________________________________________________________
######AUSTEN#######
coeff_intercept_aus <- matrix(0, m, m)
coeff_diff_aus <- matrix(0, m, m)
for (i in 1:m){
  for (j in i:m){
    coeff_intercept_aus[i, j] <- coefficients(lm(gl_proc_aus[i, j, 1:7] ~ times[17:23]))[1]
    coeff_diff_aus[i, j] <- coefficients(lm(gl_proc_aus[i, j, 1:7] ~ times[17:23]))[2]
    coeff_intercept_aus[j, i] <- coeff_intercept_aus[i, j]
    coeff_diff_aus[j, i] <- coeff_diff_aus[i, j]
  }
}

#predicted values at times in 'times_austen' are
predicted_gl_sq_B_aus <- array(0, dim=c(m,m,ntms_A))
vectorized_pred_gl_sqB_aus <- matrix(0, m * m , ntms_A)
for(i in 1:ntms_A) {
  dummy <- coeff_intercept_aus + times_austen[i] * coeff_diff_aus+rootmat(gl_proc_aus_pole)
  matched_dummy <- procOPA(rootmat(gl_proc_aus_pole), dummy, scale = FALSE, reflect = TRUE)$Bhat
  predicted_gl_sq_B_aus[, , i] <- as.matrix( proj_sparse( t(matched_dummy) %*% (matched_dummy) ) )
  vectorized_pred_gl_sqB_aus[, i] <-  as.vector(procOPA(sqrted_procmean_all,
                                                        rootmat(predicted_gl_sq_B_aus[, , i]),
                                                        scale=FALSE, reflect=TRUE)$Bhat - sqrted_procmean_all)
}


#____________________________________________________________________
######PLOT#######
plot(pca_proc$pca$x[,1:2], col=1, ylab='coordinate 2', xlab='coordinate 1', cex=0.1)
text(pca_proc$pca$x[,1:2],labels=label, col=gray.colors(25, start = 0, end = 0.7)[rank(times)], cex=0.5)
abline(lm(predict(pca_proc$pca, newdata = t(vectorized_pred_gl_sqB_dic))[,2] ~ predict(pca_proc$pca, newdata = t(vectorized_pred_gl_sqB_dic))[,1])$coef, col=1, lty=2)
abline(lm(predict(pca_proc$pca, newdata = t(vectorized_pred_gl_sqB_aus))[,2] ~ predict(pca_proc$pca, newdata = t(vectorized_pred_gl_sqB_aus))[,1])$coef, col=1)
