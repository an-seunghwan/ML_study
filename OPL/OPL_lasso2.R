rm(list=ls())
gc()
# choose seed: 3, 7, 11, 22, 30, 
set.seed(11)

if(!require("CVXR")) install.packages("CVXR")
library(CVXR)

setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/OPL")
source("OPL_lasso_init3.R")

### dimension
n = 100
p = 10
m = 5

### data setting
X = matrix(rnorm(n*p), nrow = n)
P = t(X) %*% X
P.m = matrix(apply(P, 2, mean), nrow = p, ncol = p, byrow = T)
P = P - P.m # normalize data
Aeq = matrix(runif(m*p, min = -2, max = 2), nrow = m)
beq = matrix(runif(m, min = -5, max = 5), nrow = m)

#
zero_m = matrix(rep(0, m), nrow = m)
zero_p = matrix(rep(0, p), nrow = p)
one_p = matrix(rep(1, p), nrow = p)
one_m = matrix(rep(1, m), nrow = m)
I_pp = diag(1, p)
zero_pm = matrix(rep(0, p*m), nrow = p)
zero_pp = matrix(rep(0, p*p), nrow = p)
zero_mm = matrix(rep(0, m*m), nrow = m)

### path storage
neq = dim(Aeq)[1]
maxiters = 5 * p # max number of path segments to consider
x_path = matrix(rep(0, p * maxiters), nrow = p) 
lambda_path = rep(0, maxiters) # tuning parameter
mu_patheq = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
nu_pathineq = matrix(rep(0, p * maxiters), nrow = p)
objval_path = rep(0, maxiters) # objective value

### initialization
l = init_OPL_lasso(P, Aeq, beq)
x_path[, 1] = l$x_init
lambda_path[1] = l$max_lambda
mu_patheq[, 1] = l$max_mu
nu_pathineq[, 1] = l$max_nu
active_set = l$active_set

### get direction for x and mu
H = rbind(cbind(2 * P[active_set, active_set], t(Aeq[, active_set])),
          cbind(Aeq[, active_set], zero_mm))

subgrad = -(2 * P[active_set, active_set] %*% x_path[active_set, 1] + 
              t(Aeq[, active_set]) %*% mu_patheq[, 1]) / lambda_path[1] 
# active set의 subgradient는 항상 반드시 1
# but 1이 아닌 subgradient 존재...why?

subgrad = -(2 * P %*% x_path[, 1] + t(Aeq[, ]) %*% mu_patheq[, 1]) / lambda_path[1] 
subgrad

subgrad = sign(x_path[active_set, 1, drop=F])

x_mu_dir = MASS::ginv(H) %*% rbind(subgrad, zero_m)
x_mu_dir









