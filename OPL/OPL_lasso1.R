rm(list=ls())
gc()
set.seed(11)

if(!require("CVXR")) install.packages("CVXR")
library(CVXR)

# 
n = 100
p = 10
m = 5

#
X = matrix(rnorm(n*p), nrow = n)
P = t(X) %*% X
P.m = matrix(apply(P, 2, mean), nrow = p, ncol = p, byrow = T)
P = P - P.m
Aeq = matrix(runif(m*p, min = -2, max = 2), nrow = m)
beq = matrix(runif(m, min = -5, max = 5), nrow = m)

### initialization
#
zero_m = matrix(rep(0, m), nrow = m)
zero_p = matrix(rep(0, p), nrow = p)
one_p = matrix(rep(1, p), nrow = p)
I_pp = diag(1, p)
zero_pm = matrix(rep(0, p*m), nrow = p)
zero_pp = matrix(rep(0, p*p), nrow = p)

# get initial x
# mat1 = cbind(I_pp, -I_pp)
# mat2 = cbind(Aeq, Aeq)
# mat3 = cbind(I_pp, I_pp)
# x = Variable(2 * p)
# constraints = list(mat2 %*% x == beq,
#                    mat3 %*% x >= zero_p)
# objective = Minimize(t(one_p) %*% mat1 %*% x)
# problem = Problem(objective, constraints)
# result = solve(problem)
# result$getValue(x)

# get initial x
x = Variable(p)
constraints = list(Aeq %*% x == beq,
                   x >= 0)
objective = Minimize(sum(abs(x)))
problem = Problem(objective, constraints)
result = solve(problem)
x_init = result$getValue(x)
# for numerical stability
x_init[x_init < 1e-6] = 0
x_init

# check
sum(abs(x_init))
abs(Aeq %*% x_init - beq) < 1e-6
x_init >= 0

#
# x_0 = MASS::ginv(Aeq) %*% beq

#
D1 = cbind(1, t(zero_m), t(zero_p))
D2 = cbind(zero_p, t(Aeq), I_pp)
D3 = cbind(one_p, zero_pm, zero_pp)
D4 = cbind(zero_p, zero_pm, I_pp)

#
target = Variable(1 + m + p)
constraints = list(2 * P %*% x_init + D2 %*% target <= D3 %*% target,
                   2 * P %*% x_init + D2 %*% target >= - D3 %*% target,
                   D1 %*% target >= 0,
                   D4 %*% target >= zero_p) # 제약조건 수정 필요
objective = Minimize(D1 %*% target)
problem = Problem(objective, constraints)
result = solve(problem)
result$getValue(target)
max_lambda = result$getValue(target)[1]
max_mu = result$getValue(target)[2:(1+m), , drop=F]
max_nu = result$getValue(target)[(m+2):nrow(result$getValue(target)), ,drop=F]
max_lambda; max_mu; max_nu

# boundary check
2 * P %*% x_init + t(Aeq) %*% max_mu + max_nu # every elements is almost zero






