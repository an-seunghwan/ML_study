rm(list = ls())
gc()

#####################################
# JUST FOR ORDINARY REGRESSION LOSS #
#####################################

#
setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/conlasso_eq/src")
source("conlasso_init.R")
source("conlasso_eq.R")

#
set.seed(520)

### setting --------------------------------------------------
# dimension
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

n = dim(X)[1]
p = dim(X)[2]

### equality constraints
# Aeq = matrix(rnorm(m*p, 0, 1), nrow = m)
Aeq = matrix(runif(m*p, min = -1, max = 2), nrow = m)
beq = matrix(rep(0, m), nrow = m)

# penalty weight*****(NOT USED)
penwt = rep(1, p) # 20-element Array{Frhoat64,1}

###
l = conlasso_eq(X, y, Aeq, beq, show_plotting=T)
l
