### loop part -------------------------------------------------------
k = 6
# find direction
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
if(k != 2) {
  nonzero_active_beta = which(beta_path[active_set, k-1] != 0)
  delta_rho = min(beta_path[which(active_set)[nonzero_active_beta], k-1] / delta_b[nonzero_active_beta, ])
  # predictor shrink to zero
  sub_viol = which.min(beta_path[which(active_set)[nonzero_active_beta], k-1] / delta_b[nonzero_active_beta, ])
  delta_rho_vec = c(delta_rho_vec, delta_rho)
} else {
  delta_rho_vec = c(delta_rho_vec, -1)
}

# 2. dual feasibility & KKT - stationarity condition
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

# predictor on boundary(for update active set)
delta_rho_tmp = result$getValue(delta)
kkt_viol1 = which(abs(t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta_rho_tmp * delta_b)) + 
                        t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta_rho_tmp * delta_l) - 
                        (rho_path[k-1] - delta_rho_tmp) * matrix(rep(1, sum(!active_set)), ncol = 1)) <= 1e-6)
kkt_viol2 = which(abs(t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta_rho_tmp * delta_b)) + 
                        t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta_rho_tmp * delta_l) + 
                        (rho_path[k-1] - delta_rho_tmp) * matrix(rep(1, sum(!active_set)), ncol = 1)) <= 1e-6)
kkt_viol = sort(union(kkt_viol1, kkt_viol2))

# choose delta_rho
if(delta_rho_vec[1] < 0) {
  delta_rho = delta_rho_vec[2]
  min_idx = 2
} else {
  delta_rho = min(delta_rho_vec)
  min_idx = which.min(delta_rho_vec)
}

### updates
# case 1) subgradient rule violated
if(min_idx == 1) {
  rho_path[k] = rho_path[k-1] - delta_rho
  beta_path[active_set, k] = beta_path[active_set, k-1, drop=F] - delta_rho * delta_b
  objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
  lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_l
  active_set[sub_viol] = F
  num_active = sum(active_set)
} else if(min_idx == 2) {
  rho_path[k] = rho_path[k-1] - delta_rho
  beta_path[active_set, k] = beta_path[active_set, k-1, drop=F] - delta_rho * delta_b
  objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
  lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_l
  active_set[which(!active_set)[kkt_viol]] = T
  num_active = sum(active_set)
}
