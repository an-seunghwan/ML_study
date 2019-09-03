rm(list=ls())
gc()
set.seed(520)
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")
# setting
n = 100
p = 2
m = 2

# data
X = matrix(rnorm(n*p), nrow = n)
X.m = apply(X, 2, mean)
X.m = matrix(X.m, nrow = n, ncol = p, byrow = T)
X = X - X.m
true_b = rep(1, p)
y = X %*% true_b + rnorm(n)
y = y - mean(y)

# constraints
# Aeq = matrix(sample(seq(-1, 2, by = 1), m * p, replace = T), nrow = m)
Aeq = matrix(c(1,1,
               -1,1), nrow = m, byrow = T)
beq = matrix(rep(0, m), nrow = m)

# solve
target_result = list()
for(i in 1:p) {
  Aeq1 = matrix(c(0, 0, t(Aeq[, i])), nrow = 1)
  Aeq2 = matrix(c(0, 1, rep(0, m)), nrow = 1)
  Aeq3 = matrix(c(1, rep(0, m+1)), nrow = 1)
  
  target = Variable(1 + 1 + m)
  constraints = list(Aeq2 %*% target == Aeq1 %*% target,
                     Aeq2 %*% target <= -t(X[, i]) %*% y + Aeq3 %*% target,
                     Aeq2 %*% target >= -t(X[, i]) %*% y - Aeq3 %*% target,
                     Aeq3 %*% target >= 0)
  objective = Minimize(Aeq3 %*% target)
  problem = Problem(objective, constraints)
  result = solve(problem)
  target_result[[i]] = result$getValue(target)
}
rho = sapply(target_result, function(e) {e[1, ]})
rho_max = max(sapply(target_result, function(e) {e[1, ]}))
z = sapply(target_result, function(e) {e[2, ]})
lambda = sapply(target_result, function(e) {e[3:(2+m), ]})
rho; rho_max; z; lambda

# check
diag(t(Aeq) %*% lambda)
abs(z - diag(t(Aeq) %*% lambda))
idx = which.max(rho)

-t(X) %*% y



