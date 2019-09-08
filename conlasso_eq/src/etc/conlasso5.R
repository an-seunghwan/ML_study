rm(list = ls())
gc()

################################
# JUST FOR ORDINARY REGRESSION #
################################

#
setwd("C:/Users/dpelt/OneDrive - 서울시립대학교/Documents/GitHub/ML_study/conlasso_eq/src")
source("conlasso_init.R")

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
Aeq = matrix(rnorm(m*p, 0, 1), nrow = m)
beq = matrix(rep(0, m), nrow = m)

# penalty weight
penwt = rep(1, p) # 20-element Array{Frhoat64,1}

# variable definition
neq = dim(Aeq)[1]
maxiters = 5 * p # max number of path segments to consider
beta_path = matrix(rep(0, p * maxiters), nrow = p) 
lambda_patheq = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
rho_path = rep(0, maxiters) # tuning parameter
objval_path = rep(0, maxiters) # objective value
violation_path = rep(Inf, maxiters) 

### initialization --------------------------------------------------
H = t(X) %*% X 
# find the maximum ρ (starting value)
l = path_init(X, y, Aeq, beq) # no inequality constraints
rho_path[1] = l$rho_max
lambda_patheq[, 1] = l$lambda_max

#
length(l$activeset) == m + 1

# subgradient for init_beta = 0
beta_path[, 1] = 0
resid = y - X %*% beta_path[, 1]
subgrad = (- t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]) / rho_path[1]

# active set
active_set = rep(F, p)
active_set[l$activeset] = T
num_active = sum(active_set)

### loop part -------------------------------------------------------
k = 3
# find direction
H = t(X) %*% X
M = rbind(cbind(H[active_set, active_set], t(Aeq[, active_set])),
          cbind(Aeq[, active_set], matrix(rep(0, m*m), nrow = m)))
S = rbind(subgrad[active_set, , drop=F], matrix(rep(0, m), ncol = 1))
b_m = MASS::ginv(M) %*% S
if(min(eigen(M)$values) > 0) {
  b_m = solve(M, S)
}
delta_b = b_m[1:num_active, ,drop=F]
delta_m = b_m[(num_active+1):nrow(b_m), ,drop=F]

# find delta_rho
delta_rho_vec = c()
# 1. active -> inactive
delta_rho = min(beta_path[active_set, 1] / delta_b)
delta_rho_vec = c(delta_rho_vec, delta_rho)

# 2. subgradient
#########################################3
# subgradient 조건 없이도 제대로 작동함
# 식 확인 필요
delta_rho = Inf
for(j in 1:num_active) {
  idx = which(active_set)[j]
  if(abs(subgrad[active_set, ][j] - 1) < 1e-4) {
    delta = Variable(1)
    constraints = list((-t(X[, idx]) %*% (y - X[, idx, drop=F] %*% (beta_path[idx, k-1] - delta * delta_b[j]))
                        - t(Aeq[, idx]) %*% (lambda_patheq[, k-1] - delta * delta_m)) < (rho_path[1] - delta),
                       delta >= 0)
    objective = Maximize(delta)
    problem = Problem(objective, constraints)
    result = solve(problem)
    delta_rho = min(delta_rho, result$value)
  }
  if(abs(subgrad[active_set, ][j] + 1) < 1e-4) {
    delta = Variable(1)
    constraints = list((-t(X[, idx]) %*% (y - X[, idx, drop=F] %*% (beta_path[idx, k-1] - delta * delta_b[j]))
                        - t(Aeq[, idx]) %*% (lambda_patheq[, k-1] - delta * delta_m)) > -(rho_path[1] - delta),
                       delta >= 0)
    objective = Maximize(delta)
    problem = Problem(objective, constraints)
    result = solve(problem)
    delta_rho = min(delta_rho, result$value)
  }
}
delta_rho_vec = c(delta_rho_vec, delta_rho)

# 3. dual feasibility + stationarity condition
delta = Variable(1)
constraints = list(t(Aeq[, !active_set]) %*% (lambda_patheq[, 1] - delta * delta_m) <=
                     -t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta * delta_b)) +
                     (rho_path[k-1] - delta) * matrix(rep(1, sum(!active_set)), ncol = 1),
                   t(Aeq[, !active_set]) %*% (lambda_patheq[, 1] - delta * delta_m) >=
                     -t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta * delta_b)) -
                     (rho_path[k-1] - delta) * matrix(rep(1, sum(!active_set)), ncol = 1),
                   delta >= 0)
objective = Maximize(delta)
problem = Problem(objective, constraints)
result = solve(problem)
delta_rho_vec = c(delta_rho_vec, result$value)
# predictor on boundary
delta_rho = result$value
new_j1 = which(abs(t(Aeq[, !active_set]) %*% (lambda_patheq[, 1] - delta_rho * delta_m) - 
  (-t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta_rho * delta_b)) +
  (rho_path[k-1] - delta_rho) * matrix(rep(1, sum(!active_set)), ncol = 1))) <= 1e-4)
new_j2 = which(abs(t(Aeq[, !active_set]) %*% (lambda_patheq[, 1] - delta_rho * delta_m) - 
  (-t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta_rho * delta_b)) -
  (rho_path[k-1] - delta_rho) * matrix(rep(1, sum(!active_set)), ncol = 1))) <= 1e-4)
new_j = union(new_j1, new_j2)
# pick one delta_rho
delta_rho_vec # 지나치게 작은 값은 빼면 안되나?
delta_rho = min(delta_rho_vec[which(0 < delta_rho_vec)])

# update
delta_rho = result$value
beta_path[active_set, k] = beta_path[active_set, k-1] - delta_rho * delta_b
lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_m
rho_path[k] = rho_path[k-1] - delta_rho

Aeq %*% beta_path[, k]

