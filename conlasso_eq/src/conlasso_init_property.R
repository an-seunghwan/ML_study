rm(list = ls())
gc()

################################
# JUST FOR ORDINARY REGRESSION #
################################

#
setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/conlasso_eq/src")
source("conlasso_init.R")

#
set.seed(520)

### setting --------------------------------------------------
# dimension
n = 200
p = 20
m = 10

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
Aeq = matrix(rnorm(m*p, 0, 1), nrow = m)
beq = matrix(rep(0, m), nrow = m)

# penalty weight
penwt = rep(1, p) # 20-element Array{Frhoat64,1}

# variable definition
neq = dim(Aeq)[1]
maxiters = 5 * p # max number of path segments to consider
beta_path = matrix(rep(0, p * maxiters), nrow = p) 
lambda_patheq = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
rho_path = rep(0, maxiters) # tuning parameter
objval_path = rep(0, maxiters) # objective value
violation_path = rep(Inf, maxiters) 

### initialization --------------------------------------------------
H = t(X) %*% X 
# find the maximum ρ (starting value)
l = path_init(X, y, Aeq, beq) # no inequality constraints
rho_path[1] = l$rho_max
lambda_patheq[, 1] = l$lambda_max

#
length(l$activeset) == m + 1

# subgradient for init_beta = 0
beta_path[, 1] = 0
resid = y - X %*% beta_path[, 1]
subgrad = (- t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]) / rho_path[1]

# active set
active_set = rep(F, p)
active_set[l$activeset] = T
num_active = sum(active_set)

### loop part -------------------------------------------------------
k = 2
# find direction
H = t(X) %*% X
M = rbind(cbind(H[active_set, active_set], t(Aeq[, active_set])),
          cbind(Aeq[, active_set], matrix(rep(0, m*m), nrow = m)))
S = rbind(subgrad[active_set, , drop=F], matrix(rep(0, m), ncol = 1))
b_m = MASS::ginv(M) %*% S
if(min(eigen(M)$values) > 0) {
  b_m = solve(M, S)
}
delta_b = b_m[1:num_active, ,drop=F]
delta_m = b_m[(num_active+1):nrow(b_m), ,drop=F]

# 1) whole active set
set = l$activeset
Z = solve(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set])) # invertible
Matrix::rankMatrix(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set])) == 
  ncol(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set]))
W = t(Aeq[, set]) %*% Z %*% Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) # NOT identity matrix
W
abs(sum(diag(W)) - m) < 1e-6
length(set) - sum(diag(W))
abs(length(set) - sum(diag(W)) - 1) < 1e-6

# 2) b == m
set = l$activeset[1:m]
Z = solve(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set])) # invertible
Matrix::rankMatrix(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set])) == 
  ncol(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set]))
Matrix::rankMatrix(Aeq[, set]) == length(set) # column is linearly independent
W = t(Aeq[, set]) %*% Z %*% Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) # identity matrix
W
length(set) - sum(diag(W))

# 3) b < m
set = l$activeset[1:3]
Z = solve(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set])) # NOT invertible
Z = MASS::ginv(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set]))
Matrix::rankMatrix(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set])) == 
  ncol(Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) %*% t(Aeq[, set]))
Matrix::rankMatrix(Aeq[, set]) == length(set) # column is linearly independent
W = t(Aeq[, set]) %*% Z %*% Aeq[, set] %*% solve(t(X[, set]) %*% X[, set]) # identity matrix
W
length(set) - sum(diag(W))
