rm(list = ls())
gc()

#
setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/classo")
source("classopath_init.R")
source("constrsparsereg.R")
source("classopath_all.R")
source("classopath_one.R")
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")

#############################################################################
# Calculate the solution path of the constrained lasso problem that minimizes
# `0.5sumabs2(√obswt .* (y - X * β)) + ρ * sumabs(penwt .* β)`
# subject to linear constraints.
#############################################################################

set.seed(520)

# setting
n = 500
p = 20
m = 10

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

#
# X = data.matrix(longley[, 1:6])
# n = dim(X)[1]
# p = dim(X)[2]
# y = data.matrix(longley[, 7])
# X.m = apply(X, 2, mean)
# X.m = matrix(X.m, nrow = n, ncol = p, byrow = T)
# X = X - X.m
# y = y - mean(y)
# m = 3

# n = dim(X)[1]
# p = dim(X)[2]

### equality constraints
# default
# Aeq = matrix(0, nrow = 0, ncol = dim(X)[2])
# beq = rep(0, dim(Aeq)[1])
# use
Aeq = matrix(sample(seq(-2, 2, by = 1), m * p, replace = T), nrow = m)
beq = matrix(rep(0, m), nrow = m)

### inequality constraints
# default
Aineq = matrix(0, nrow = 0, ncol = dim(X)[2])
bineq = rep(0, dim(Aineq)[1])
# use 
# Aineq = -diag(rep(1, p))
# bineq = rep(0, p)

# penalty weight
penwt = rep(1, p) # 20-element Array{Frhoat64,1}

# 1. init 단계에서 violation된 모든 beta를 activate 하는 경우
result_all = classopath_all(X, y, Aeq, beq, Aineq, bineq, penwt)

# result
cat("Test: ", abs(Aeq %*% result_all$beta_path[, result_all$steps]), "\n")
beta_all = result_all$beta_path[, result_all$steps]

# plot
beta_max = max(result_all$beta_path)
beta_min = min(result_all$beta_path)
sum_abs_beta = apply(abs(result_all$beta_path), 2, sum)
#
# par(mfrow = c(1,1))
# plot(sum_abs_beta, result_all$beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
#      xlab = "sum(abs(beta))", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) case one')
# for(i in 2:p) points(sum_abs_beta, result_all$beta_path[i, ], type = 'l', col = i, lwd = 2)
#
par(mfrow = c(1,1))
plot(seq(1, result_all$steps), result_all$beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
     xlab = "steps", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) case one')
for(i in 2:p) points(seq(1, result_all$steps), result_all$beta_path[i, ], type = 'l', col = i, lwd = 2)

# 2. init 단계에서 초기 잔차와 가장 correlation이 높은 predictor 1개만을 선택
result_one = classopath_one(X, y, Aeq, beq, Aineq, bineq, penwt)

# result
cat("Test: ", abs(Aeq %*% result_one$beta_path[, result_one$steps]), "\n")
beta_one = result_one$beta_path[, result_one$steps]

# plot
beta_max = max(result_one$beta_path)
beta_min = min(result_one$beta_path)
sum_abs_beta = apply(abs(result_one$beta_path), 2, sum)
#
# par(mfrow = c(1,1))
# plot(sum_abs_beta, result_one$beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
#      xlab = "sum(abs(beta))", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) case two')
# for(i in 2:p) points(sum_abs_beta, result_one$beta_path[i, ], type = 'l', col = i, lwd = 2)
#
par(mfrow = c(1,1))
plot(seq(1, result_one$steps), result_one$beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
     xlab = "steps", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) case two')
for(i in 2:p) points(seq(1, result_one$steps), result_one$beta_path[i, ], type = 'l', col = i, lwd = 2)

#
par(mfrow = c(1,1))
min = min(beta_all, beta_one)
max = max(beta_all, beta_one)
max_step = max(result_all$steps, result_one$steps)
plot(seq(1, max_step), ylim = c(min, max), type = 'l', col = 1, lwd = 2,
     xlab = "steps", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) compare')
for(i in 1:p) points(seq(1, result_all$steps), result_all$beta_path[i, ], type = 'l', col = 1, lwd = 2)
for(i in 1:p) points(seq(1, result_one$steps), result_one$beta_path[i, ], type = 'l', col = 2, lwd = 2, lty = 2)

# compare
plot(beta_all - beta_one, ylim = c(-1, 1), type = "h")
beta_all
beta_one
# exactly same solution
# KKT condition, subgradient rule, sign mismatch 검사만 제대로 하면 맨처음 active set은 상관 없음