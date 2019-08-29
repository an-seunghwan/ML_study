# rm(list=ls())
# gc()
# set.seed(520)
# if(!require("CVXR")) install.packages("CVXR")
# library("CVXR")
# # setting
# n = 200
# p = 10
# m = 5
# 
# # data
# X = matrix(rnorm(n*p), nrow = n)
# X.m = apply(X, 2, mean)
# X.m = matrix(X.m, nrow = n, ncol = p, byrow = T)
# X = X - X.m
# true_b = rep(1, p)
# y = X %*% true_b + rnorm(n)
# y = y - mean(y)
# 
# # constraints
# Aeq = matrix(sample(seq(-2, 2, by = 1), m * p, replace = T), nrow = m)
# beq = matrix(rep(0, m), nrow = m)

# 
init_path = function(X, y, Aeq, beq) {
  # beq should be all zero vector
  n = dim(X)[1]
  p = dim(X)[2]
  m = dim(Aeq)[1]
  
  # setting
  zero_p = matrix(rep(0, p), ncol = 1)
  zero_m = matrix(rep(0, m), ncol = 1)
  zero_pm = matrix(rep(0, p*m), nrow = p)
  zero_pp = matrix(rep(0, p*p), nrow = p)
  one_p = matrix(rep(1, p), ncol = 1)
  one_m = matrix(rep(1, m), ncol = 1)
  I_p = matrix(diag(rep(1, p)), nrow = p)
  
  # design matrix for Constraints
  D1 = matrix(c(1, t(zero_p), t(zero_m)), nrow = 1)
  D2 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(zero_p, I_p, zero_pm),
             c(0, t(zero_p), t(zero_m)))
  D3 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(zero_p, zero_pp, t(Aeq)),
             c(0, t(zero_p), t(zero_m)))
  D4 = rbind(0, -t(X) %*% y, 0)
  D5 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(one_p, zero_pp, zero_pm),
             c(0, t(zero_p), t(zero_m)))
  
  # solve
  target = Variable(1 + p + m)
  constraints = list(D2 %*% target == D3 %*% target,
                     D2 %*% target <= D4 + D5 %*% target,
                     D2 %*% target >= D4 - D5 %*% target,
                     D1 %*% target >= 0)
  objective = Minimize(D1 %*% target)
  problem = Problem(objective, constraints)
  result = solve(problem)
  
  # result
  target = result$getValue(target)
  rho_max = target[1, ]
  # z = target[2:(1+p), ,drop=F]
  lambda_max = target[(p+2):nrow(target), ,drop=F]
  
  # violation
  idx1 = which(abs((-t(X) %*% y + rho_max * one_p) - t(Aeq) %*% lambda_max) <= 1e-4)
  idx2 = which(abs((-t(X) %*% y - rho_max * one_p) - t(Aeq) %*% lambda_max) <= 1e-4)
  activeset = sort(union(idx1, idx2))
  
  # full rank check
  flag = Matrix::rankMatrix(t(Aeq)[activeset, ]) == ncol(t(Aeq)[activeset, ])
  
  return(list(rho_max = rho_max,
              lambda_max = lambda_max,
              activeset = activeset,
              flag = flag))
}

# l = init_path(X, y, Aeq, beq)
# l$rho_max
# l$lambda_max
# l$activeset
# l$flag