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

#
one_p = matrix(rep(1, p), nrow = p)
rho = 0
lambda_min = MASS::ginv(t(Aeq)) %*% (-t(X) %*% y - rho * one_p)
lambda_max = MASS::ginv(t(Aeq)) %*% (-t(X) %*% y + rho * one_p)

lambda = (lambda_max + lambda_min) / 2
t(Aeq) %*% lambda
(-t(X) %*% y - rho * one_p <= t(Aeq) %*% lambda) * (-t(X) %*% y + rho * one_p >= t(Aeq) %*% lambda) 

# 
rho = 100000000
while(T) {
  lambda_min = MASS::ginv(t(Aeq)) %*% (-t(X) %*% y - rho * one_p)
  lambda_max = MASS::ginv(t(Aeq)) %*% (-t(X) %*% y + rho * one_p)
  
  min_flag = (-t(X) %*% y - rho * one_p <= t(Aeq) %*% lambda_min)
  max_flag = (-t(X) %*% y + rho * one_p >= t(Aeq) %*% lambda_max)
  
  if(sum(min_flag & max_flag) == p) break
  rho = rho + 1
}
rho

# solve
rho = 1000
target = Variable(1 + m)
Amat1 = matrix(c(rep(0, m), 1), nrow = 1)
Amat = rbind(cbind(t(Aeq), one_p), cbind(-t(Aeq), one_p))
bmat = rbind(-t(X) %*% y + 2 * rho * one_p, t(X) %*% y + 2 * rho * one_p)
constraints = list(Amat %*% target <= bmat, Amat1 %*% target <= rho)
objective = Maximize(Amat1 %*% target)
problem = Problem(objective, constraints)
result = solve(problem)
lambda = result$getValue(target)[1:m, , drop=F]
lambda
z = result$getValue(target)[m+1, ]
z

#
rho = 0 # rho is increasing direction
threshold1 = 0.1
threshold2 = 0.1
while(T) {
  target = Variable(1 + m)
  Amat1 = matrix(c(rep(0, m), 1), nrow = 1)
  Amat = rbind(cbind(t(Aeq), one_p), cbind(-t(Aeq), one_p))
  bmat = rbind(-t(X) %*% y + 2 * rho * one_p, t(X) %*% y + 2 * rho * one_p)
  constraints = list(Amat %*% target <= bmat, Amat1 %*% target <= rho)
  objective = Maximize(Amat1 %*% target)
  problem = Problem(objective, constraints)
  result = solve(problem)
  z = result$getValue(target)[m+1, ]
  lambda = result$getValue(target)[1:m, , drop=F]
  # cat("z: ", z, "rho: ", rho, "\n")
  
  violation_check = (-t(X) %*% y - rho * one_p <= t(Aeq) %*% lambda) * (-t(X) %*% y + rho * one_p >= t(Aeq) %*% lambda)
  if(sum(violation_check) == p) break
  # if(abs(z - rho) < threshold1) break
  rho_old = rho
  lambda_old = lambda
  z_old = z
  rho = rho + threshold2
}
violation_check = (-t(X) %*% y - rho_old * one_p <= t(Aeq) %*% lambda_old) * (-t(X) %*% y + rho_old * one_p >= t(Aeq) %*% lambda_old)
violation_idx = which(violation_check == 0)
lambda
z_old
rho
violation_idx

# [1,] -57.304791
# [2,] -40.777656
# [3,] -11.697450
# [4,]  21.200860
# [5,]  -6.103915

#
lambda = result$getValue(target)[1:m, , drop=F]
lambda
z
rho
violation_check
violation_idx
(-t(X) %*% y - rho * one_p <= t(Aeq) %*% lambda) * (-t(X) %*% y + rho * one_p >= t(Aeq) %*% lambda)



(-t(X) %*% y - rho * one_p <= t(Aeq) %*% lambda) * (-t(X) %*% y + rho * one_p >= t(Aeq) %*% lambda) 
t(Aeq) %*% lambda

# result
target = result$getValue(target)
rho_min = target[1, ]
z = target[2:(1+p), ,drop=F]
lambda = target[(p+2):nrow(target), ,drop=F]
rho_min; z; lambda
Aeq

# check
sum(abs(z - t(Aeq) %*% lambda)) < 1e-6

# violation check : rho is decreasing
t(Aeq) %*% lambda
((-t(X) %*% y + rho_min) - (t(Aeq) %*% lambda)) < 1e-8
which.min((-t(X) %*% y + rho_min) - (t(Aeq) %*% lambda))

# maximum of lowerbound
-t(X) %*% y + rho_min
-t(X) %*% y - rho_min
min(-t(X) %*% y + rho_min)
max(-t(X) %*% y - rho_min)

t(Aeq) %*% lambda

t(Aeq[, 1]) %*% lambda