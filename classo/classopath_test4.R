rm(list = ls())
gc()

#
setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/classo")
source("classopath_init.R")
source("classopath3.R")
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")

#############################################################################
# Calculate the solution path of the constrained lasso problem that minimizes
# `0.5sumabs2(√obswt .* (y - X * β)) + ρ * sumabs(penwt .* β)`
# subject to linear constraints.
#############################################################################

set.seed(520)

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

n = dim(X)[1]
p = dim(X)[2]

### equality constraints
# default
# Aeq = matrix(0, nrow = 0, ncol = dim(X)[2])
# beq = rep(0, dim(Aeq)[1])
# use
# Aeq = matrix(sample(seq(-1, 1, by = 1), m * p, replace = T), nrow = m)
Aeq = matrix(rnorm(m*p, 1, 1), nrow = m)
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

##################################################################################
# 1. init 단계에서 violation된 모든 beta를 activate 하는 경우
result_all = classopath_modified(X, y, Aeq, beq, Aineq, bineq, penwt, choose_one = F)

# result
cat("Test: ", abs(Aeq %*% result_all$beta_path[, result_all$steps]), "\n")
beta_all = result_all$beta_path[, result_all$steps]

# plot
# par(mfrow = c(1,1))
# plot(sum_abs_beta, result_all$beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
#      xlab = "sum(abs(beta))", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) case one')
# for(i in 2:p) points(apply(abs(result_all$beta_path), 2, sum), result_all$beta_path[i, ], type = 'l', col = i, lwd = 2)

par(mfrow = c(1,1))
plot(seq(1, result_all$steps), result_all$beta_path[1, ], 
     ylim = c(min(result_all$beta_path), max(result_all$beta_path)), type = 'l', col = 1, lwd = 2,
     xlab = "steps", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) all')
for(i in 2:p) points(seq(1, result_all$steps), result_all$beta_path[i, ], type = 'l', col = i, lwd = 2)

###############################################################################
# 2. init 단계에서 초기 잔차와 가장 correlation이 높은 predictor 1개만을 선택
result_one = classopath_modified(X, y, Aeq, beq, Aineq, bineq, penwt, choose_one = T)

# result
cat("Test: ", abs(Aeq %*% result_one$beta_path[, result_one$steps]), "\n")
beta_one = result_one$beta_path[, result_one$steps]

# plot
# par(mfrow = c(1,1))
# plot(sum_abs_beta, result_one$beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
#      xlab = "sum(abs(beta))", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) case two')
# for(i in 2:p) points(apply(abs(result_one$beta_path), 2, sum), result_one$beta_path[i, ], type = 'l', col = i, lwd = 2)

par(mfrow = c(1,1))
plot(seq(1, result_one$steps), result_one$beta_path[1, ], 
     ylim = c(min(result_one$beta_path), max(result_one$beta_path)), type = 'l', col = 1, lwd = 2,
     xlab = "steps", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) one')
for(i in 2:p) points(seq(1, result_one$steps), result_one$beta_path[i, ], type = 'l', col = i, lwd = 2)

############################################################################
# compare
par(mfrow = c(1,1))
min = min(result_all$beta_path, result_one$beta_path)
max = max(result_all$beta_path, result_one$beta_path)
max_step = max(result_all$steps, result_one$steps)
plot(seq(1, max_step), ylim = c(min, max), type = 'l', col = 1, lwd = 2,
     xlab = "steps", ylab = "beta", main = 'BETA coefs PATH(Constrained LASSO) compare')
for(i in 1:p) points(seq(1, result_all$steps), result_all$beta_path[i, ], type = 'l', col = 1, lwd = 2)
for(i in 1:p) points(seq(1, result_one$steps), result_one$beta_path[i, ], type = 'l', col = 2, lwd = 2, lty = 2)

#
plot(beta_all - beta_one, ylim = c(-0.1, 0.1), type = "h")
beta_all - beta_one
beta_all
beta_one
# exactly same solution

#
t(Aeq)
sum(abs(Aeq %*% result_all$beta_path[, result_all$steps])) # all
sum(abs(Aeq %*% result_one$beta_path[, result_one$steps])) # one
result_all$beta_path
result_one$beta_path
##################################################################################################
# KKT condition, subgradient rule, sign mismatch 검사만 제대로 하면 맨처음 active set은 상관 없음#
# 왜냐면 첫번째 업데이트 이전에 검사하게 되면 active set이 동일해짐
# 한번씩 필요한 step가 다른 것 빼고는 동일
# BUT 제약조건을 만족시키는 정도에서 차이남, 근데 완전 같을때도 있음
##################################################################################################