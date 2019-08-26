rm(list = ls())
gc()

#
setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/classo")
source("find_rho_max.R")
source("constrsparsereg.R")
source("classopath.R")
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")

#############################################################################
# Calculate the solution path of the constrained lasso problem that minimizes
# `0.5sumabs2(√obswt .* (y - X * β)) + ρ * sumabs(penwt .* β)`
# subject to linear constraints.
#############################################################################

### truth with sum constrant sum(beta) = 0
# β = zeros(p)
# β[1:round(Int, p / 4)] = 0
# β[(round(Int, p / 4) + 1):round(Int, p / 2)] = 1
# β[(round(Int, p / 2) + 1):round(Int, 3p / 4)] = 0
# β[(round(Int, 3p / 4) + 1):p] = -1

### generate data - normalize된 data 사용
X=as.matrix(read.csv("C:/Julia/classo/src/X1.csv", header = F)) # n * p
y=as.matrix(read.csv("C:/Julia/classo/src/y1.csv", header = F)) # n

n = dim(X)[1]
p = dim(X)[2]

### equality constraints
# default
# Aeq = matrix(0, nrow = 0, ncol = dim(X)[2])
# beq = rep(0, dim(Aeq)[1])
# use
Aeq = matrix(rep(1, p), nrow = 1) # 1×20 Array{Frhoat64,2}
beq = matrix(0, nrow = 1) # 1-element Array{Fsrhoat64,1}
# another example
# Aeq = matrix(1, nrow = 2, p); Aeq[1,11:20] = 0; Aeq[2,1:10] = 0
# beq = matrix(0, nrow = 2, 1)

### inequality constraints
# default
Aineq = matrix(0, nrow = 0, ncol = dim(X)[2])
bineq = rep(0, dim(Aineq)[1])
# use 
# Aineq = -diag(rep(1, p))
# bineq = rep(0, p)

# penalty weight
penwt = rep(1, p) # 20-element Array{Frhoat64,1}

# do it!
result = classopath(X, y, Aeq, beq, Aineq, bineq, penwt)

# result
cat("Test: sum of beta < 1e-6 =", sum(result$beta_path[, result$steps]) < 1e-6, "\n")

# plot
par(mfrow = c(1,1))
beta_max = max(result$beta_path)
beta_min = min(result$beta_path)
sum_abs_beta = apply(abs(result$beta_path), 2, sum)
plot(sum_abs_beta, result$beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
     xlab = "sum(abs(beta))", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO)')
for(i in 2:p) points(sum_abs_beta, result$beta_path[i, ], type = 'l', col = i, lwd = 2)
