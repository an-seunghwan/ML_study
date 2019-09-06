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
# Aeq = matrix(sample(c(-2,-1,0,1,2), m*p, replace = T), nrow = m)
beq = matrix(rep(0, m), nrow = m)

# penalty weight*****
penwt = rep(1, p) # 20-element Array{Frhoat64,1}

# variable definition
neq = dim(Aeq)[1]
maxiters = 5 * p # max number of path segments to consider
beta_path = matrix(rep(0, p * maxiters), nrow = p) 
lambda_patheq = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
rho_path = rep(0, maxiters) # tuning parameter
objval_path = rep(0, maxiters) # objective value
# violation_path = rep(Inf, maxiters) 

### initialization --------------------------------------------------
H = t(X) %*% X 
# find the maximum ρ (starting value)
l = path_init(X, y, Aeq, beq) # no inequality constraints
rho_path[1] = l$rho_max
lambda_patheq[, 1] = l$lambda_max

# subgradient for init_beta = 0
beta_path[, 1] = 0
resid = y - X %*% beta_path[, 1]
subgrad = (- t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]) / rho_path[1]
objval_path[1] = sum((y - X %*% beta_path[ ,1])^2) + rho_path[1] * sum(abs(beta_path[1]))

# active set
active_set = rep(F, p)
active_set[l$activeset] = T
num_active = sum(active_set)

### main path following algorithm -------------------------------------------
for(k in 2:maxiters) {
  ### find direction
  H = t(X) %*% X
  M = rbind(cbind(H[active_set, active_set], t(Aeq[, active_set])),
            cbind(Aeq[, active_set], matrix(rep(0, m*m), nrow = m)))
  S = rbind(subgrad[active_set, , drop=F], matrix(rep(0, m), ncol = 1))
  b_l = MASS::ginv(M) %*% S
  if(min(eigen(M)$values) > 0) {
    b_l = solve(M, S)
  }
  delta_b = b_l[1:num_active, ,drop=F]
  delta_l = b_l[(num_active+1):nrow(b_l), ,drop=F]
  
  ### find delta_rho
  delta_rho_vec = c()
  # 1. active -> inactive (== subgradient rule check)
  if(k != 2) { # at initial, all beta is zero
    nonzero_active_beta = which(beta_path[active_set, k-1] != 0)
    delta_rho = min(beta_path[which(active_set)[nonzero_active_beta], k-1] / delta_b[nonzero_active_beta, ])
    # predictor shrink to zero
    sub_viol = which.min(beta_path[which(active_set)[nonzero_active_beta], k-1] / delta_b[nonzero_active_beta, ])
    delta_rho_vec = c(delta_rho_vec, delta_rho)
    
  } else {
    delta_rho_vec = c(delta_rho_vec, -1)
  }
  
  # 2. dual feasibility & KKT - stationarity condition
  if(sum(active_set) != p) { # if all predictors are activated, do not run opt problem
    delta = Variable(1)
    constraints = list(t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta * delta_b)) + 
                         t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta * delta_l) 
                       <= (rho_path[k-1] - delta) * matrix(rep(1, sum(!active_set)), ncol = 1),
                       
                       t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta * delta_b)) +
                         t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta * delta_l) 
                       >= -(rho_path[k-1] - delta) * matrix(rep(1, sum(!active_set)), ncol = 1),
                       
                       delta >= 0)
    objective = Maximize(delta)
    problem = Problem(objective, constraints)
    result = solve(problem)
    delta_rho_vec = c(delta_rho_vec, result$getValue(delta))
    
    # predictor on boundary(for update active set) -> candidates for active_set
    delta_rho_tmp = result$getValue(delta)
    kkt_viol1 = which(abs(t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta_rho_tmp * delta_b)) + 
                            t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta_rho_tmp * delta_l) - 
                            (rho_path[k-1] - delta_rho_tmp) * matrix(rep(1, sum(!active_set)), ncol = 1)) <= 1e-6)
    kkt_viol2 = which(abs(t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta_rho_tmp * delta_b)) + 
                            t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta_rho_tmp * delta_l) + 
                            (rho_path[k-1] - delta_rho_tmp) * matrix(rep(1, sum(!active_set)), ncol = 1)) <= 1e-6)
    kkt_viol = sort(union(kkt_viol1, kkt_viol2))
    
  } else {
    delta_rho_vec = c(delta_rho_vec, -1)
  }
  
  # choose delta_rho
    # 1. KKT condition is only violated
  if(delta_rho_vec[1] < 0 & delta_rho_vec[2] > 0) {
    delta_rho = delta_rho_vec[2]
    min_idx = 2
    # 2. subgradient rule is only violated
  } else if(delta_rho_vec[1] > 0 & delta_rho_vec[2] < 0){
    delta_rho = delta_rho_vec[1]
    min_idx = 1
    # 3. choose first violated rule
  } else if(delta_rho_vec[1] > 0 & delta_rho_vec[2] > 0){
    delta_rho = min(delta_rho_vec)
    min_idx = which.min(delta_rho_vec)
    # 4. there is no candidates for delta_rho -> termination of algorithm
  } else if(delta_rho_vec[1] < 0 & delta_rho_vec[2] < 0){
    delta_rho = rho_path[k-1] # make rho to zero
    min_idx = 0
  }
  
  ### updates
    # case 1) subgradient rule is violated
  if(min_idx == 1) {
    rho_path[k] = rho_path[k-1] - delta_rho
    beta_path[active_set, k] = beta_path[active_set, k-1, drop=F] - delta_rho * delta_b
    objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
    lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_l
    active_set[sub_viol] = F
    num_active = sum(active_set)
    # case 2) KKT condition is violated
  } else if(min_idx == 2) {
    rho_path[k] = rho_path[k-1] - delta_rho
    beta_path[active_set, k] = beta_path[active_set, k-1, drop=F] - delta_rho * delta_b
    objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
    lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_l
    active_set[which(!active_set)[kkt_viol]] = T
    num_active = sum(active_set)
    # case 3) terminate algorithm
  } else if(min_idx == 0) {
    rho_path[k] = rho_path[k-1] - delta_rho
    beta_path[active_set, k] = beta_path[active_set, k-1, drop=F] - delta_rho * delta_b
    objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
    lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_l
    break
  }
}

# result
rho_path[1:k]
beta_path[, 1:k]
lambda_patheq[, 1:k]
objval_path[1:k]

# equality constraint check
sum(abs(Aeq %*% beta_path[, k])) <= 1e-6
