library(CorporaCoCo)
library(Matrix)

label <- c("TTC", "BR", "BH", "DC", "DS", "GE", "HT", "LD", "MC", "NN", "OT",
           "OMF", "PP", "ED", "OCS", "C", 
           "EM", "PE", "PR", "LS", "MA", "NO", "SE" )

author_col <- c(rep(1,16), rep(2,7))

times <- c(1859, 1841, 1852, 1849, 1846, 1860, 1854, 1855, 1843, 1838, 1837,
           1864, 1836, 1870, 1840, 1843,
           1814, 1815, 1796, 1794, 1811, 1798, 1795)

#loading novels in
novel_list <- list()
for (i in 1:23){
  novel_list[[i]] <- readRDS(paste0('austen_and_dickens', '/', label[i], '.rds'))
}

#wordlist
lookup <- sort( table( unlist(novel_list) ), decreasing=TRUE )
lookup <- rownames(lookup)
words <- lookup[1:1000]


#creating adjacency matrices of co-occurences
# you must use sparse matrices else you will run out of memory
novel_adj <- lapply(novel_list, function(x) {
  s <- surface(x, span = "5LR")  # from CorporaCoCo (see help)
  sparseMatrix(
    i = match(s$x, lookup),
    j = match(s$y, lookup),
    x = s$H,
    dims = rep(length(lookup), 2)
  )
})
names(novel_adj) <- c()


#creating graph Laplacians
List_LC_full <- list() #weighted gl
List_L_full <- list() #weighted gl scaled by trace
List_LC <- list() #weighted gl truncated for memory
List_L <- list() #weighted gl scaled by trace truncated for memory
L <- array(0, dim=c(1000, 1000, 23)) #list L as an array
LC <- array(0, dim=c(1000, 1000, 23))# list LC as an array
for (i in 1:23){
  List_LC_full[[i]] <- -(novel_adj [i])[[1]]
  diag(List_LC_full[[i]] ) <- 0
  diag(List_LC_full[[i]]) <- -rowSums(List_LC_full[[i]])
  List_L_full[[i]] <- List_LC_full[[i]]/sum(diag(List_LC_full[[i]]))
  #
  List_LC[[i]] <- as.matrix(List_LC_full[[i]][1:1000,1:1000])
  diag(List_LC[[i]]) <- 0
  diag(List_LC[[i]]) <- -rowSums(List_LC[[i]])
  List_L[[i]] <- List_LC[[i]]/sum(diag(List_LC[[i]]))
  #
  LC[,,i] <- List_LC[[i]]
  L[,,i] <- List_L[[i]]
}

remove(i, novel_list, novel_adj)
