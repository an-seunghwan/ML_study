rm(list=ls())
gc()
set.seed(520)
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")
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

# constraints
Aeq = matrix(sample(seq(-1, 2, by = 1), m * p, replace = T), nrow = m)
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
Aeq

abs((-t(X) %*% y + rho_min * one_p) - (t(Aeq) %*% lambda)) # right boundary
which.min(abs((-t(X) %*% y + rho_min * one_p) - (t(Aeq) %*% lambda)))

# find violation
threshold = 1e-8
rho = rho_min
count = 1
while (T) {
  rho = rho - threshold
  violation_check = (-t(X) %*% y - rho * one_p <= t(Aeq) %*% lambda) * (-t(X) %*% y + rho * one_p >= t(Aeq) %*% lambda)
  violation_idx = which(violation_check == 0)
  if(sum(violation_check) < p) break
  count = count + 2
}
count
violation_idx
lambda
rho_min

#
# (-t(X) %*% y - rho * one_p)[violation_idx, ] 
# (t(Aeq) %*% lambda)[violation_idx, ] 
# (-t(X) %*% y + rho * one_p)[violation_idx, ] 

#
abs((-t(X) %*% y + rho_min * one_p) - (t(Aeq) %*% lambda))







