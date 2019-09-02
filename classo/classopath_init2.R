rm(list = ls())
gc()

#
setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/classo")
source("classopath_init.R")
source("classopath3.R")
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")
# if(!require("cape")) install.packages("cape")
# library("cape")

#############################################################################
# Calculate the solution path of the constrained lasso problem that minimizes
# `0.5sumabs2(√obswt .* (y - X * β)) + ρ * sumabs(penwt .* β)`
# subject to linear constraints.
#############################################################################

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
Aeq = matrix(sample(c(1,2,0), m * p, replace = T), nrow = m)
# Aeq = matrix(rnorm(m*p, 1, 1), nrow = m)
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

# subgradient for init_beta = 0
beta_path[, 1] = 0
resid = y - X %*% beta_path[, 1]
subgrad = - t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]

# choose candidates
# findDepMat(Aeq[, l$activeset])
Aeq[, l$activeset]
l$activeset
setB = c(1,2,3,5,9)
setB = l$activeset

# Matrix::rankMatrix(Aeq[, l$activeset])
Aeq[, setB]
Matrix::rankMatrix(Aeq[, setB]) == length(setB)

new_setActive = rep(F, p)
new_setActive[setB] = T
# derivative for beta and lambda
M = rbind(cbind(t(X[, new_setActive]) %*% X[, new_setActive], t(Aeq[, new_setActive])),
          cbind(Aeq[, new_setActive], matrix(rep(0, neq * neq), nrow = neq)))
delta_beta_lambda = MASS::ginv(M) %*% rbind(subgrad[new_setActive, ,drop=F] / rho_path[1], 
                                            matrix(rep(0, neq), ncol = 1))
delta_beta = delta_beta_lambda[1:sum(new_setActive), ,drop=F]
delta_lambda = delta_beta_lambda[(sum(new_setActive)+1):nrow(delta_beta_lambda), ,drop=F]

###############################################################################################
t(X[, new_setActive]) %*% X[, new_setActive]
t(Aeq[, new_setActive]) %*% solve(diag(rep(.Machine$double.eps, neq))) %*% Aeq[, new_setActive]
t(X[, new_setActive]) %*% X[, new_setActive] - 
  t(Aeq[, new_setActive]) %*% solve(diag(rep(.Machine$double.eps, neq))) %*% Aeq[, new_setActive]

MASS::ginv(t(X[, new_setActive]) %*% X[, new_setActive] - 
        t(Aeq[, new_setActive]) %*% solve(diag(rep(.Machine$double.eps, neq))) %*% Aeq[, new_setActive])
MASS::ginv(t(X[, new_setActive]) %*% X[, new_setActive] - 
             t(Aeq[, new_setActive]) %*% solve(diag(rep(.Machine$double.eps, neq))) %*% Aeq[, new_setActive]) %*% 
  subgrad[new_setActive, ,drop=F] / rho_path[1]

MASS::ginv(t(X[, new_setActive]) %*% X[, new_setActive])
MASS::ginv(t(X[, new_setActive]) %*% X[, new_setActive]) %*% subgrad[new_setActive, ,drop=F] / rho_path[1]

# check constraints
abs(t(X[, !new_setActive]) %*% X[, new_setActive] %*% delta_beta +
      t(Aeq[, !new_setActive]) %*% delta_lambda) <= matrix(rep(1, sum(!new_setActive)), ncol = 1)
abs(t(X[, !new_setActive]) %*% y +
    t(Aeq[, !new_setActive]) %*% lambda_patheq[, 1]) <= matrix(rep(rho_path[1], sum(!new_setActive)), ncol = 1)

#####################################################################################
if(!l$flag) return(0) # check lambda is unique
setActive = matrix(rep(F, p), ncol = 1)
setActive[l$activeset, ] = T
nActive = sum(setActive != 0)

beta_path[, 1] = 0
objval_path[1] = sum(t(X) %*% y) / 2 # init_beta is zero

#
mu_pathineq[mu_pathineq < 0] = 0
residIneq = Aineq %*% beta_path[, 1] - bineq
setIneqBorder = residIneq == 0
nIneqBorder = sum(setIneqBorder != 0)

# initialize subgradient vector (stationarity condition)
resid = y - X %*% beta_path[, 1]
if(neq > 0 & nineq > 0) {
  subgrad = - t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1] - t(Aineq) %*% mu_pathineq[, 1]
} else if(neq > 0 & nineq == 0) {
  subgrad = - t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]
} else if(neq == 0 & nineq > 0) {
  subgrad = - t(X) %*% resid - t(Aineq) %*% mu_pathineq[, 1]
}

# subgradient of beta
# subgrad[setActive] = sign(beta_path[setActive, 1])
# subgrad[!setActive] = subgrad[!setActive] / rho_path[1]
subgrad = subgrad / rho_path[1] # 처음에는 beta가 다 0이니까 sign말고 직접 계산?

# calculate degrees of freedom
# rankAeq = Matrix::rankMatrix(Aeq) ### Brian's comment: need to make it more efficient
rankAeq = ifelse(dim(Aeq)[1] > 0, Matrix::rankMatrix(Aeq), 0)
df_path[1] = nActive - rankAeq - nIneqBorder

# set initial violations counter to 0
violation_path[1] = 0

# sign for path direction (originally went both ways, but increasing was retired -> only decreasing)
dirsgn = -1




















