rm(list = ls())
gc()

if(!require("CVXR")) install.packages("CVXR")
library("CVXR")

#####################################
# JUST FOR ORDINARY REGRESSION LOSS #
#####################################

#
set.seed(520)

### setting --------------------------------------------------
# dimension
n = 100
p = 10
k = 5

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

### inequality constraints
# Aeq = matrix(rnorm(m*p, 0, 1), nrow = m)
Aineq = matrix(runif(k*p, min = -1, max = 2), nrow = k)
bineq = matrix(rep(0, k), nrow = k)

# penalty weight*****(NOT USED)
penwt = rep(1, p) # 20-element Array{Frhoat64,1}

n = dim(X)[1]
p = dim(X)[2]
k = dim(Aineq)[1]

# setting
zero_p = matrix(rep(0, p), ncol = 1)
zero_k = matrix(rep(0, k), ncol = 1)
zero_pk = matrix(rep(0, p*k), nrow = p)
zero_kp = matrix(rep(0, p*k), nrow = k)
zero_pp = matrix(rep(0, p*p), nrow = p)
one_p = matrix(rep(1, p), ncol = 1)
one_k = matrix(rep(1, k), ncol = 1)
I_p = matrix(diag(rep(1, p)), nrow = p)
I_k = matrix(diag(rep(1, k)), nrow = k)

# design matrix for Constraints
D1 = matrix(c(1, t(zero_p), t(zero_k)), nrow = 1)
D2 = rbind(c(0, t(zero_p), t(zero_k)),
           cbind(zero_p, I_p, zero_pk),
           c(0, t(zero_p), t(zero_k)))
D3 = rbind(c(0, t(zero_p), t(zero_k)),
           cbind(zero_p, zero_pp, t(Aineq)),
           c(0, t(zero_p), t(zero_k)))
D4 = rbind(0, t(X) %*% y, 0)
D5 = rbind(c(0, t(zero_p), t(zero_k)),
           cbind(one_p, zero_pp, zero_pk),
           c(0, t(zero_p), t(zero_k)))
D6 = cbind(zero_k, zero_kp, I_k)

# solve
target = Variable(1 + p + k)
constraints = list(D2 %*% target == D3 %*% target,
                   D2 %*% target <= D4 + D5 %*% target,
                   D2 %*% target >= D4 - D5 %*% target,
                   D1 %*% target >= 0,
                   D6 %*% target >= zero_k)
objective = Minimize(D1 %*% target)
problem = Problem(objective, constraints)
result = solve(problem)

# result
target = result$getValue(target)
rho_max = target[1, ]
# z = target[2:(1+p), ,drop=F]
nu_max = target[(p+2):nrow(target), ,drop=F]

# violation
idx1 = which(abs(-t(X) %*% y + rho_max * one_p + t(Aineq) %*% nu_max) <= 1e-4)
idx2 = which(abs(-t(X) %*% y - rho_max * one_p + t(Aineq) %*% nu_max) <= 1e-4)
activeset = sort(union(idx1, idx2))

# full rank check
flag = T
if(dim(Aineq)[1] > 1) flag = Matrix::rankMatrix(t(Aineq)[activeset, ]) == ncol(t(Aineq)[activeset, ])