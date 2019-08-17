rm(list=ls())
gc()
library(mvtnorm)

X = read.table("C:/Users/dpelt/OneDrive/바탕 화면/Mayson/UOS_graduate/Convex optimiaztion/gmm_data.txt", sep = "\t")
X = cbind(X$V1, X$V2)

gmmcluster = function(X, num_cluster, plotting = F) {
  # init parameter
  init = function()
  {
    l = list()
    l = lapply(1:num_cluster, function(i) {
      list(weight = 1 / num_cluster,
           means = runif(2),
           cov = matrix(c(1, 0, 0, 1), nrow=2))
    })
    return(l)
  }
  l = init()
  
  n = nrow(X)
  max_iter = 100
  
  for(iter in 1:max_iter)
  {
    ### 1. E step ###
    # responsibility
    r = sapply(l, function(e) {
      e$weight * dmvnorm(X, e$means, e$cov)
    })
    sum_r = apply(r, 1, sum)
    e_r = r / matrix(sum_r, nrow = n, ncol = num_cluster, byrow = F)
    num_k = apply(e_r, 2, sum)
    
    ### 2. M step ###
    # mean
    mu = lapply(1:num_cluster, function(i) {
      t(e_r[, i]) %*% X / num_k[i]
    })
    # cov
    sigma = lapply(1:num_cluster, function(i) {
      temp_mu = X - matrix(mu[[i]], nrow = n, ncol = ncol(X), byrow = T)
      (t(matrix(e_r[ , i], nrow = n, ncol = ncol(X), byrow = F) * temp_mu) %*% temp_mu) / num_k[i]
    })
    # pi
    pi = num_k / n
    # update parameters
    l = lapply(1:num_cluster, function(i) {
      list(weight = pi[i], means = mu[[i]], cov = sigma[[i]])
    })
  }
  
  weights = sapply(l, function(e) { e$weight })
  means = sapply(l, function(e) { e$means })
  covs = sapply(l, function(e) { e$cov })

  if(plotting) {
    plot(X)
    points(t(means), col=2:(num_cluster + 1), cex=2)
    xpts = seq(range(X[,1])[1], range(X[,1])[2], length.out=50)
    ypts = seq(range(X[,2])[1], range(X[,2])[2], length.out=50)
    contour_grid = as.matrix(expand.grid(xpts,ypts))
    for(i in 1:num_cluster)
    {
      mix_contour = weights[i] * dmvnorm(contour_grid, t(means[,i]), matrix(covs[,i], nrow = 2, byrow = T))
      contour(xpts, ypts, matrix(mix_contour,length(xpts),length(ypts)), add = T, col = i+1)
    }
  }
  return(list(weight = sapply(l, function(e) { e$weight }),
              means = sapply(l, function(e) { e$means }),
              cov = sapply(l, function(e) { e$cov })))
}

parameter = gmmcluster(X, num_cluster = 6, plotting = T)  
