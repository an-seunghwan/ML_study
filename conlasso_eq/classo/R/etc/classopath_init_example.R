rm(list=ls())
gc()
set.seed(520)
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")
# setting
n = 100
p = 5
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
Aeq = matrix(sample(c(-2, -1, 1, 2), m * p, replace = T), nrow = m)
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

idx1 = which(abs((-t(X) %*% y + rho_min * one_p) - t(Aeq) %*% lambda) <= 1e-4)
idx2 = which(abs((-t(X) %*% y - rho_min * one_p) - t(Aeq) %*% lambda) <= 1e-4)
activeset = union(idx1, idx2)
sort(activeset)

# plot
l1 = seq(-100, 200, length.out = 1000)
l2 = seq(-100, 200, length.out = 1000)

# rho_min = rho_min - 100

y_plus = function(l1, i) (-t(X) %*% y + rho_min * one_p)[i, ] / t(Aeq)[i ,2] - (l1 * t(Aeq)[i ,1]) / t(Aeq)[i ,2]
y_minus = function(l1, i) (-t(X) %*% y - rho_min * one_p)[i, ] / t(Aeq)[i ,2] - (l1 * t(Aeq)[i ,1]) / t(Aeq)[i ,2]

plot(l1, l2, ylab = "", type = "n", col = 1)
points(lambda[1, ], lambda[2, ], col = p+1)
for(i in 1:p) lines(l1, y_plus(l1, i), ylab = "", col = i)
for(i in 1:p) lines(l1, y_minus(l1, i), ylab = "", col = i)
# 1: 검은색, 3: 초록색, 5: 청록색 (겹치는 직선 색깔)