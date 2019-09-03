rm(list = ls())
gc()

#
setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/classo")
source("classopath_init.R")
source("classopath3.R")
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")

#
set.seed(520)

# setting
n = 100
p = 10
m = 5

# data
X = matrix(rnorm(n*p), nrow = n)
X.m = apply(X, 2, mean)
X.m = matrix(X.m, nrow = n, ncol = p, byrow = T)
X = X - X.m
true_b = rep(1, p)
y = X %*% true_b + rnorm(n)
y = y - mean(y)

n = dim(X)[1]
p = dim(X)[2]

### equality constraints
# default
# Aeq = matrix(0, nrow = 0, ncol = dim(X)[2])
# beq = rep(0, dim(Aeq)[1])
# use
# Aeq = matrix(sample(seq(-1, 1, by = 1), m * p, replace = T), nrow = m)
# Aeq = matrix(sample(c(1,2,0), m * p, replace = T), nrow = m)
Aeq = matrix(rnorm(m*p, 0, 1), nrow = m)
# Aeq = matrix(c(1,1,1,0,0,0,1,1,1,0,
# -1,-1,-1,0,0,0,-1,-1,-1,0,
# 0,0,0,0,0,2,2,2,2,2), nrow = m, byrow = T)
# 0,0,0,0,0,-1,-1,-1,-1,-1,
# 1,1,1,1,1,1,1,1,1,1,
# 0,0,1,2,1,2,1,0,0,2,
# 0,0,2,4,2,4,2,0,0,4), nrow = m, byrow = T)
beq = matrix(rep(0, m), nrow = m)

### inequality constraints
# default
Aineq = matrix(0, nrow = 0, ncol = dim(X)[2])
bineq = rep(0, dim(Aineq)[1])
# use 
# Aineq = -diag(rep(1, p))
# bineq = rep(0, p)

# penalty weight
penwt = rep(1, p) # 20-element Array{Frhoat64,1}
penidx = rep(T, p)
choose_one = F

# alrhocate variables arhong path
neq = dim(Aeq)[1]
nineq = dim(Aineq)[1]
maxiters = 5 * (p + nineq) # max number of path segments to consider
beta_path = matrix(rep(0, p * maxiters), nrow = p) 
lambda_patheq = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
mu_pathineq = matrix(rep(0, nineq * maxiters), nrow = nineq) # dual variables for inequality
rho_path = rep(0, maxiters) # tuning parameter
df_path = rep(Inf, maxiters) # degree of freedom
objval_path = rep(0, maxiters) # objective value
violation_path = rep(Inf, maxiters) 

### initialization
H = t(X) %*% X 
# find the maximum ρ (starting value)
l = init_path(X, y, Aeq, beq) # no inequality constraints
rho_path[1] = l$rho_max
lambda_patheq[, 1] = l$lambda_max

#
length(l$activeset) == m + 1

# subgradient for init_beta = 0
beta_path[, 1] = 0
resid = y - X %*% beta_path[, 1]
subgrad = (- t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]) / rho_path[1]

#
active_set = rep(F, p)
active_set[l$activeset] = T
# active_set[c(2,3,5,7,9,13,15)] = T
# active_set[5] = T

Matrix::rankMatrix(t(X[, active_set]) %*% X[, active_set]) == ncol(t(X[, active_set]) %*% X[, active_set])
XX_inv = solve(t(X[, active_set]) %*% X[, active_set])

Matrix::rankMatrix(Aeq[, active_set])

Aeq[, active_set] %*% XX_inv %*% t(Aeq[, active_set])
Matrix::rankMatrix(Aeq[, active_set] %*% XX_inv %*% t(Aeq[, active_set])) == ncol(Aeq[, active_set] %*% XX_inv %*% t(Aeq[, active_set]))
# AXXA_inv = solve(Aeq[, active_set] %*% XX_inv %*% t(Aeq[, active_set]))
solve(Aeq[, active_set] %*% XX_inv %*% t(Aeq[, active_set]))
AXXA_inv = solve(Aeq[, active_set] %*% XX_inv %*% t(Aeq[, active_set]))

AXXA_inv = MASS::ginv(Aeq[, active_set] %*% XX_inv %*% t(Aeq[, active_set]))

t(Aeq[, active_set]) %*% AXXA_inv %*% Aeq[, active_set]
t(X[, active_set]) %*% X[, active_set]

Z = t(Aeq[, active_set]) %*% AXXA_inv %*% Aeq[, active_set] %*% XX_inv
Z
diag(Z)
















