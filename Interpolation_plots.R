source('functions_needed.R')
load("workspace for novels.RData")

x <- List_L
m<-1000


####____________________________________________________________________


#Euclidean
dickens_eucmean <- Mean_GL(x[1:16], euc=TRUE, sqrt=FALSE, proc=FALSE)
austen_eucmean <- Mean_GL(x[17:23], euc=TRUE, sqrt=FALSE, proc=FALSE)
int_euc_0pt5 <- Interpolation(dickens_eucmean, austen_eucmean,0.5, euc=TRUE, sqrt=FALSE, proc=FALSE)$int
int_euc_6 <-Interpolation(dickens_eucmean, austen_eucmean, 6, euc=TRUE, sqrt=FALSE, proc=FALSE)$int
int_euc_min5 <-Interpolation(dickens_eucmean, austen_eucmean, -5, euc=TRUE, sqrt=FALSE, proc=FALSE)$int

rownames(int_euc_0pt5)<-lookup[1:m]
colnames(int_euc_0pt5)<-lookup[1:m]
rownames(int_euc_6)<-lookup[1:m]
colnames(int_euc_6)<-lookup[1:m]
rownames(int_euc_min5)<-lookup[1:m]
colnames(int_euc_min5)<-lookup[1:m]


ggraph_plot(int_euc_0pt5[1:25,1:25],0, quantile_edges=0.95,  tol=-1*10^-16, colour_node = 'black', colour='black')
ggraph_plot(int_euc_6[1:25,1:25],0, quantile_edges=0.95,  tol=-1*10^-16, colour_node = 'black', colour='black')
ggraph_plot(int_euc_min5[1:25,1:25],0, quantile_edges=0.95,  tol=-1*10^-16, colour_node = 'black', colour='black')


####__________________________________________________________________


#square root Euclidean
dickens_sqrtmean <- proj_sparse( Mean_GL(x[1:16], euc=FALSE, sqrt=TRUE, proc=FALSE) )
austen_sqrtmean  <- proj_sparse( Mean_GL(x[17:23], euc=FALSE, sqrt=TRUE, proc=FALSE) )

#int_sqrt_min9 <-Interpolation(dickens_sqrtmean, austen_sqrtmean, -9, euc=FALSE, sqrt=TRUE, proc=FALSE)$int
int_sqrt_0pt5 <-Interpolation(dickens_sqrtmean, austen_sqrtmean, 0.5, euc=FALSE, sqrt=TRUE, proc=FALSE)$int
int_sqrt_6 <-Interpolation(dickens_sqrtmean, austen_sqrtmean, 6, euc=FALSE, sqrt=TRUE, proc=FALSE)$int
int_sqrt_min5 <-Interpolation(dickens_sqrtmean, austen_sqrtmean, -5, euc=FALSE, sqrt=TRUE, proc=FALSE)$int

rownames(int_sqrt_0pt5)<-lookup[1:m]
colnames(int_sqrt_0pt5)<-lookup[1:m]
rownames(int_sqrt_6)<-lookup[1:m]
colnames(int_sqrt_6)<-lookup[1:m]
rownames(int_sqrt_min5)<-lookup[1:m]
colnames(int_sqrt_min5)<-lookup[1:m]


ggraph_plot(int_sqrt_0pt5[1:25,1:25],0, 0.95,  tol=-1*10^-16, colour_node = 'black', colour_line='black')
ggraph_plot(int_sqrt_6[1:25,1:25],0, 0.95,  tol=-1*10^-16, colour_node = 'black', colour_line='black')
ggraph_plot(int_sqrt_min5[1:25,1:25],0, 0.95,  tol=-1*10^-16, colour_node = 'black', colour_line='black')

####____________________________________________________________________


#procrustes
dickens_procmean <- proj_sparse( Mean_GL(x[1:16], euc=FALSE, sqrt=FALSE, proc=TRUE) )
austen_procmean  <- proj_sparse( Mean_GL(x[17:23], euc=FALSE, sqrt=FALSE, proc=TRUE) )

int_proc_0pt5 <-Interpolation(dickens_procmean, austen_procmean, 0.5, euc=FALSE, sqrt=FALSE, proc=TRUE)$int
int_proc_6 <-Interpolation(dickens_procmean, austen_procmean, 6, euc=FALSE, sqrt=FALSE, proc=TRUE)$int
int_proc_min5 <-Interpolation(dickens_procmean, austen_procmean, -5, euc=FALSE, sqrt=FALSE, proc=TRUE)$int


rownames(int_proc_0pt5)<-lookup[1:m]
colnames(int_proc_0pt5)<-lookup[1:m]
rownames(int_proc_6)<-lookup[1:m]
colnames(int_proc_6)<-lookup[1:m]
rownames(int_proc_min5)<-lookup[1:m]
colnames(int_proc_min5)<-lookup[1:m]


ggraph_plot(int_proc_0pt5[1:25,1:25],0, 0.95,  tol=-1*10^-16, colour_node = 'black', colour_line='black')
ggraph_plot(int_proc_min5[1:25,1:25],0, 0.95,  tol=-1*10^-16, colour_node = 'black', colour_line='black')
ggraph_plot(int_proc_6[1:25,1:25],0, 0.95,  tol=-1*10^-16, colour_node = 'black', colour_line='black')
