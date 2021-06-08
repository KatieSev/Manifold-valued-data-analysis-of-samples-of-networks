#source('functions_needed.R')
#load("workspace for novels.RData")

m<-dim(List_L[[1]])[1]
C <- 100 #permutations  (used 1000 in the paper)

########Euclidean###########
dickens_austen_euc_2 <- (dist_sq_means(L[,,1:16],L[,,17:23],euc=TRUE, sqrt=FALSE, proc=FALSE))

dist_vec_euc <- rep(0,C)
set.seed(14)
for (j in 1:C){
  rndm_samp1 <- sample(1:23, 16)
  rndm_samp2 <- setdiff(1:23, rndm_samp1)
  dist_vec_euc[j] <-dist_sq_means(L[, , rndm_samp1], L[, , rndm_samp2],
                                  euc = TRUE, sqrt = FALSE, proc = FALSE)
}

plot(density(dist_vec_euc), main='', xlab="", xlim=c(0,0.003),lwd=2) #, xlim=c(0,0.0007))
abline(v=dickens_austen_euc_2, col=2,lwd=2)
mtext(expression({d[E]^{"2"}}(bar(A), bar(B))), side=1, line=1.2, cex=0.7)


########Square root Euclidean###########
dickens_austen_sqrt_2 <- dist_sq_means(L[,,1:16],L[,,17:23],euc=FALSE, sqrt=TRUE, proc=FALSE)
dist_vec_sqrt <- rep(0,C)   
set.seed(14)
for (j in 1:C)
{
  rndm_samp1 <- sample(1:23, 16)
  rndm_samp2 <- setdiff(1:23, rndm_samp1)
  dist_vec_sqrt[j] <- dist_sq_means(L[, , rndm_samp1], L[, , rndm_samp2],
                                    euc = FALSE, sqrt = TRUE, proc = FALSE)
}

plot(density(dist_vec_sqrt), main='', xlab=expression({d[H]^{"2"}} (bar(A), bar(B))), xlim=c(0,0.03))
abline(v=dickens_austen_sqrt_2, col=2)


########Procrustes###########
dickens_austen_proc_2 <- dist_sq_means(L[,,1:16],L[,,17:23],euc=FALSE, sqrt=FALSE, proc=TRUE)
dist_vec_proc <- rep(0,C)  
set.seed(14)
for (j in 1:C)
{
  rndm_samp1 <- sample(1:23, 16)
  rndm_samp2 <- setdiff(1:23, rndm_samp1)
  dist_vec_proc[j] <- dist_sq_means(L[, , rndm_samp1], L[, , rndm_samp2],
                                    euc = FALSE, sqrt = FALSE, proc = TRUE)
}
plot(density(dist_vec_proc), main='', xlab=expression({d[S]^{"2"}} (bar(A), bar(B))), xlim=c(0,0.025))
abline(v=dickens_austen_proc_2,col=2)
