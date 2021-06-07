#source('functions_needed.R')
#load("workspace for novels.RData")

#####################
#dickens
euc_mean <-Mean_GL(L[,,1:16])

sqrt_mean <-Mean_GL(L[,,1:16],euc=FALSE, sqrt = TRUE)

proc_mean <-Mean_GL(L[,,1:16], euc=FALSE, proc=TRUE)


#####################
#austen

euc_mean <-Mean_GL(L[,,17:23])

sqrt_mean <-Mean_GL(L[,,17:23],euc=FALSE, sqrt = TRUE)

proc_mean <-Mean_GL(L[,,17:23], euc=FALSE, proc=TRUE)
