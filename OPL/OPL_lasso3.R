rm(list=ls())
gc()
# choose seed: 3, 7, 11, 22, 30, 
set.seed(11)

if(!require("CVXR")) install.packages("CVXR")
library(CVXR)

setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/OPL")
source("OPL_lasso_init3.R")

### dimension
n = 100
k = 4
p = 3*k
m = 5

#
zero_m = matrix(rep(0, m), nrow = m)
zero_p = matrix(rep(0, p), nrow = p)
one_p = matrix(rep(1, p), nrow = p)
one_m = matrix(rep(1, m), nrow = m)
I_pp = diag(1, p)
zero_pm = matrix(rep(0, p*m), nrow = p)
zero_pp = matrix(rep(0, p*p), nrow = p)
zero_mm = matrix(rep(0, m*m), nrow = m)

### data setting
X = matrix(rnorm(n*p, sd=0.1), nrow = n)
P = t(X) %*% X
P.m = matrix(apply(P, 2, mean), nrow = p, ncol = p, byrow = T)
P = P - P.m # normalize data
# mu_t = matrix(runif(p, min=0.5, max=2), ncol = 1)
# q_t = matrix(runif(p, min=1, max=2), ncol = 1)
# rho_t = runif(1, min=3, max=4)
# w_hat = matrix(rnorm(p, mean = 0.5, sd=0.2), ncol = 1)
# Aeq = rbind(cbind(t(mu_t), t(q_t), -t(q_t)),
#             cbind(t(one_p), t(zero_p), t(zero_p)),
#             cbind(I_pp, -I_pp, I_pp))
# beq = matrix(c(rho_t, 1, w_hat), ncol = 1)
a = runif(k)
Aeq = rbind(cbind(t(runif(k, max=2)), t(a), -t(a)),
            cbind(t(rep(1, k)), t(rep(0, k*2))),
            cbind(diag(k), -diag(k), diag(k)))
w_hat = c(0.2, 0.5, 0.1, 0.2)
beq = matrix(c(runif(1, min=0, max=1), 1, w_hat), ncol = 1)

### path storage
neq = dim(Aeq)[1]
maxiters = 5 * p # max number of path segments to consider
x_path = matrix(rep(0, p * maxiters), nrow = p) 
lambda_path = rep(0, maxiters) # tuning parameter
mu_patheq = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
nu_pathineq = matrix(rep(0, p * maxiters), nrow = p)
objval_path = rep(0, maxiters) # objective value

### initialization
l = init_OPL_lasso(P, Aeq, beq)
x_path[, 1] = l$x_init
lambda_path[1] = l$max_lambda
mu_patheq[, 1] = l$max_mu
nu_pathineq[, 1] = l$max_nu
active_set = l$active_set
num_active = sum(active_set)

### get direction for x and mu
H = rbind(cbind(2 * P[active_set, active_set], t(Aeq[, active_set])),
          cbind(Aeq[, active_set], zero_mm))

subgrad = -(2 * P[active_set, active_set] %*% x_path[active_set, 1] + 
              t(Aeq[, active_set]) %*% mu_patheq[, 1]) / lambda_path[1] 
subgrad
# active set의 subgradient는 항상 반드시 1
# but 1이 아닌 subgradient 존재...why?
# 심지어 -1도 있음

-(2 * P %*% x_path[, 1, drop=F] + t(Aeq) %*% mu_patheq[, 1, drop=F]) / lambda_path[1]
# subgrad

# subgrad = sign(x_path[active_set, 1, drop=F])
# subgrad

x_mu_dir = MASS::ginv(H) %*% rbind(subgrad, zero_m)
x_mu_dir
delta_x = x_mu_dir[1:num_active, , drop=F]
delta_mu = x_mu_dir[(num_active+1):nrow(x_mu_dir), , drop=F]
delta_x; delta_mu

Aeq[, active_set]
# t(Aeq[, active_set]) %*% x_mu_dir[(sum(active_set)+1):nrow(x_mu_dir), , drop=F]
abs(t(Aeq[ , active_set]) %*% delta_mu - subgrad) < 1e-6








