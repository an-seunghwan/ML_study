rm(list=ls())
gc()
set.seed(520)
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")
# setting
n = 200
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

# constraints
Aeq = matrix(sample(seq(-1, 2, by = 1), m * p, replace = T), nrow = m)
# Aeq = matrix(sample(c(-2, -1, 1, 2), m * p, replace = T), nrow = m)
# Aeq = matrix(c(1,1,0,
# -1,0,1,
# -1,-1,1), nrow = m, byrow = T)
beq = matrix(rep(0, m), nrow = m)

# setting
zero_p = matrix(rep(0, p), ncol = 1)
zero_m = matrix(rep(0, m), ncol = 1)
zero_pm = matrix(rep(0, p*m), nrow = p)
zero_pp = matrix(rep(0, p*p), nrow = p)
one_p = matrix(rep(1, p), ncol = 1)
one_m = matrix(rep(1, m), ncol = 1)
I_p = matrix(diag(rep(1, p)), nrow = p)

# design matrix for Constraints
Aeq1 = matrix(c(1, t(zero_p), t(zero_m)), nrow = 1)
Aeq2 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(zero_p, I_p, zero_pm),
             c(0, t(zero_p), t(zero_m)))
Aeq3 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(zero_p, zero_pp, t(Aeq)),
             c(0, t(zero_p), t(zero_m)))
Aeq4 = rbind(0, -t(X)%*%y, 0)
Aeq5 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(one_p, zero_pp, zero_pm),
             c(0, t(zero_p), t(zero_m)))

# solve
target = Variable(1 + p + m)
constraints = list(Aeq2 %*% target == Aeq3 %*% target,
                   Aeq2 %*% target <= Aeq4 + Aeq5 %*% target,
                   Aeq2 %*% target >= Aeq4 - Aeq5 %*% target,
                   Aeq1 %*% target >= 0)
objective = Minimize(Aeq1 %*% target)
problem = Problem(objective, constraints)
result = solve(problem)

# result
target = result$getValue(target)
rho_min = target[1, ]
z = target[2:(1+p), ,drop=F]
lambda = target[(p+2):nrow(target), ,drop=F]
rho_min; z; lambda
t(Aeq)

# rho_min에서는 모든 부등조건을 만족
idx1 = which((-t(X) %*% y + rho_min * one_p) <= t(Aeq) %*% lambda)
idx2 = which((-t(X) %*% y - rho_min * one_p) >= t(Aeq) %*% lambda)
activeset = c(idx1, idx2)
sort(activeset)

### 1. 이건 경계에 있는 predictor를 바로 activeset 지정 -> 이건 아닌듯!
idx1 = which(abs((-t(X) %*% y + rho_min * one_p) - t(Aeq) %*% lambda) <= 1e-4)
idx2 = which(abs((-t(X) %*% y - rho_min * one_p) - t(Aeq) %*% lambda) <= 1e-4)
activeset = union(idx1, idx2)
sort(activeset)

### 2. violation check
# find new_lambda corresponding to decreased rho(아주 작게 감소된 rho에 대해 violation을 새로 체크)
new_rho = rho_min - 1e-4 # rho is decreasing direction
target = Variable(1 + m)
Amat1 = matrix(c(rep(0, m), 1), nrow = 1)
Amat = rbind(cbind(t(Aeq), one_p), cbind(-t(Aeq), one_p))
bmat = rbind(-t(X) %*% y + 2 * new_rho * one_p, t(X) %*% y + 2 * new_rho * one_p)
constraints = list(Amat %*% target <= bmat, 
                   Amat1 %*% target <= new_rho)
objective = Maximize(Amat1 %*% target)
problem = Problem(objective, constraints)
result = solve(problem)
z = result$getValue(target)[m+1, ]
new_lambda = result$getValue(target)[1:m, , drop=F]
# 새로운 rho, lambda에 대해 부등조건을 만족하지 않는(범위를 벗어나는) predictor를 찾는다
idx1 = which((-t(X) %*% y + new_rho * one_p) <= t(Aeq) %*% new_lambda)
idx2 = which((-t(X) %*% y - new_rho * one_p) >= t(Aeq) %*% new_lambda)
activeset = c(idx1, idx2)
sort(activeset) # 최종 active set

# 3. linear combination
activeset = c()
for(i in 1:m) {
  idx1 = which(abs(((-t(X) %*% y + rho_min * one_p) - t(Aeq)[, -i, drop=F] %*% 
                      lambda[-i, ,drop=F]) / t(Aeq)[, i] - lambda[i, ]) < 1e-4)
  idx2 = which(abs(((-t(X) %*% y - rho_min * one_p) - t(Aeq)[, -i, drop=F] %*% 
                      lambda[-i, ,drop=F]) / t(Aeq)[, i] - lambda[i, ]) < 1e-4)
  # print(c(idx1, idx2))
  # if(length(activeset) > 0) activeset = c(idx1, idx2)
  # activeset = intersect(activeset, c(idx1, idx2))
  activeset = union(activeset, c(idx1, idx2))
}
sort(activeset)