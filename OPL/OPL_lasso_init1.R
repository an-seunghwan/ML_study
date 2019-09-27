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
##### solve
# min   sum(abs(beta))
# s.t.  Ax = b
#       x >= 0
#####
x = Variable(p)
constraints = list(Aeq %*% x == beq,
                   x >= 0)
objective = Minimize(sum(abs(x)))
problem = Problem(objective, constraints)
result = solve(problem)
x_init = result$getValue(x)
x_init[x_init < 1e-6] = 0.0 # for numerical stability
x_init

# conditions check
sum(abs(x_init))
abs(Aeq %*% x_init - beq) < 1e-6
x_init >= 0

#
# x_0 = MASS::ginv(Aeq) %*% beq

# Design matrix
D1 = cbind(1, t(zero_m), t(zero_p))
D2 = cbind(zero_p, t(Aeq), I_pp)
D3 = cbind(one_p, zero_pm, zero_pp)
D4 = cbind(zero_p, zero_pm, I_pp)

# get initial lambda, mu, nu corresponding to x_init
target = Variable(1 + m + p)
constraints = list(2 * P %*% x_init + D2 %*% target <= D3 %*% target, # stationarity
                   2 * P %*% x_init + D2 %*% target >= - D3 %*% target, # stationarity
                   D1 %*% target >= 0, # lambda must be positive
                   D4 %*% target >= zero_p, # dual feasibility
                   t(D4 %*% target) %*% x_init == 0 # complementary slackness
                   )
objective = Minimize(D1 %*% target)
problem = Problem(objective, constraints)
result = solve(problem)
# result$getValue(target)
max_lambda = result$getValue(target)[1]
max_mu = result$getValue(target)[2:(1+m), , drop=F]
max_nu = result$getValue(target)[(m+2):nrow(result$getValue(target)), , drop=F]
max_nu[max_nu < 1e-6] = 0.0 # for numerical stability 
x_init; max_lambda; max_mu; max_nu

# boundary check for active set
idx1 = which(abs(2 * P %*% x_init + t(Aeq) %*% max_mu + max_nu - max_lambda) < 1e-6) # Right ineq part
idx2 = which(abs(2 * P %*% x_init + t(Aeq) %*% max_mu + max_nu + max_lambda) < 1e-6) # Left ineq part
# 1. lambda를 infinite부터 감소시키는 과정에서 boundary에 놓여있는 predictor들을 선택한 acitve set
sort(union(idx1, idx2))
# 2. 초기 x값이 0이 아닌 값들을 active set으로 선택하는 방법
  # -> 이게 맞는 듯...?
which(x_init != 0)



