source('PCA_Euclidean')
m<-1000
pca<-PCA_euc$pca
cat<-author_col

times_dickens<- seq(1836,1870,1)
times_austen <- seq(1794, 1815,1)
length_tms_dickens <- length(times_dickens)
length_tms_austen <- length(times_austen)

######DICKENS#######
#Approach B
coeff_intercept_dic<-matrix(0,m,m)
coeff_diff_dic<-matrix(0,m,m)
for (i in 1:m)
{ for (j in i:m)
{coeff_dic <- coefficients(lm(L[i,j,1:16] ~ times[1:16]))
coeff_intercept_dic[i,j] <- coeff_dic[1]
coeff_diff_dic[i,j] <- coeff_dic[2]
coeff_intercept_dic[j,i] <- coeff_intercept_dic[i,j]
coeff_diff_dic[j,i] <- coeff_diff_dic[i,j]
}
}
#predited values are
predicted_glB_dic <-array(0, dim=c(m,m,length_tms_dickens))
for(i in 1:length_tms_dickens){
  predicted_glB_dic[,,i] <- as.matrix(proj_sparse( coeff_intercept_dic+times_dickens[i]*coeff_diff_dic))
  print(i)
}

vectorized_pred_glB_dic <-matrix(0,m *m ,length_tms_dickens)
for (i in 1:length_tms_dickens)
{
  vectorized_pred_glB_dic[,i] <-as.vector(predicted_glB_dic[,,i])}

plot(pca$x[,1:2], col=cat, ylab='coordinate 2', xlab='coordinate 1', cex=0.1)
text(pca$x[,1:2],labels=label, col=cat, cex=0.5)
points(predict(pca, newdata = t(vectorized_pred_glB_dic))[,1:2])
#____________________________________________________________________

########AUSTEN############
#Approach B
coeff_intercept_aus<-matrix(0,m,m)
coeff_diff_aus<-matrix(0,m,m)
for (i in 1:m)
{ for (j in i:m)
{coeff_intercept_aus[i,j] <-coefficients(lm(L[i,j,17:23] ~ times[17:23]))[1]
coeff_diff_aus[i,j] <-coefficients(lm(L[i,j,17:23] ~ times[17:23]))[2]
coeff_intercept_aus[j,i] <- coeff_intercept_aus[i,j]
coeff_diff_aus[j,i] <- coeff_diff_aus[i,j]
}
  print(i)
}
#predited values are
predicted_glB_aus <-array(0, dim=c(m,m,length_tms_austen))
for(i in 1:length_tms_austen)
  predicted_glB_aus[,,i] <- as.matrix(proj_sparse( (coeff_intercept_aus+times_austen[i]*coeff_diff_aus)))

vectorized_pred_glB_aus <-matrix(0,m *m ,length_tms_austen)
for (i in 1:length_tms_austen)
{for (j in 1:m)
  vectorized_pred_glB_aus[((j-1)*m +1):(j*m ),i] <-predicted_glB_aus[,j,i]}

plot(pca$x[,1:2], col=cat, ylab='coordinate 2', xlab='coordinate 1', cex=0.1)
text(pca$x[,1:2],labels=label, col=cat, cex=0.5)
points(predict(pca, newdata = t(vectorized_pred_glB_aus))[,1:2])

###combine####
#can also see code plotting regression on PC plot for these
plot(pca$x[,1:2], col=cat, ylab='coordinate 2', xlab='coordinate 1', cex=0.1)
text(pca$x[,1:2],labels=label, col=rainbow(25)[rank(times)], cex=0.5)
points(predict(pca, newdata = t(vectorized_pred_glB_dic))[,1:2],col=4,pch=19,cex=0.1)
points(predict(pca, newdata = t(vectorized_pred_glB_aus))[,1:2],col=2,pch=19,cex=0.1)
abline(lm(predict(pca, newdata = t(vectorized_pred_glB_dic))[,2] ~ predict(pca, newdata = t(vectorized_pred_glB_dic))[,1])$coef, col=4)
abline(lm(predict(pca, newdata = t(vectorized_pred_glB_aus))[,2] ~ predict(pca, newdata = t(vectorized_pred_glB_aus))[,1])$coef, col=2)

#or grayscale

plot(pca$x[,1:2], col=1, ylab='coordinate 2', xlab='coordinate 1', cex=0.1)
text(pca$x[,1:2],labels=label, col=gray.colors(25, start = 0, end = 0.7)[rank(times)], cex=0.5)
#points(predict(pca, newdata = t(vectorized_pred_glB_dic))[,1:2],col=1,pch=19,cex=0.1)
#points(predict(pca, newdata = t(vectorized_pred_glB_aus))[,1:2],col=1,pch=19,cex=0.1)
abline(lm(predict(pca, newdata = t(vectorized_pred_glB_dic))[,2] ~ predict(pca, newdata = t(vectorized_pred_glB_dic))[,1])$coef, col=1, lty=2)
abline(lm(predict(pca, newdata = t(vectorized_pred_glB_aus))[,2] ~ predict(pca, newdata = t(vectorized_pred_glB_aus))[,1])$coef, col=1)
