#
conlasso_eq <- function(X, y, Aeq, beq,
                        penwt = rep(1, p),
                        show_plotting = F) 
{
  n = dim(X)[1]
  p = dim(X)[2]
  
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
  # find the maximum ρ and corresponding lambda (starting value)
  l = path_init(X, y, Aeq, beq) # no inequality constraints
  rho_path[1] = l$rho_max
  lambda_patheq[, 1] = l$lambda_max
  
  # subgradient for init_beta = 0
  beta_path[, 1] = 0
  resid = y - X %*% beta_path[, 1]
  subgrad = (t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]) / rho_path[1]
  
  # loss value
  objval_path[1] = sum((y - X %*% beta_path[ ,1])^2) + rho_path[1] * sum(abs(beta_path[1]))
  
  # initial active set
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
    
    if(min(eigen(M)$values) > 0) {
      b_l = solve(M, S)
    } else {
      b_l = MASS::ginv(M) %*% S  
    }
    
    # derivative for beta and lambda for rho
    delta_b = b_l[1:num_active, ,drop=F]
    delta_l = b_l[(num_active+1):nrow(b_l), ,drop=F]
    
    ### find delta_rho
    delta_rho_vec = c()
    
    # 1. active -> inactive (== subgradient rule check)
    nonzero_active_beta = which(active_set)[which(beta_path[active_set, k-1] != 0)]
    # predictor have potential to shrink 0
      # -> which beta coef and delta_beta have different direction
    sub_mismatch = sign(beta_path[nonzero_active_beta, k-1]) * sign(delta_b[which(beta_path[active_set, k-1] != 0), ])
    
    # there exists potential
    if(any(sub_mismatch == 1)) {
      delta_rho_sub_tmp = beta_path[nonzero_active_beta[which(sub_mismatch == 1)], k-1] / 
        delta_b[which(beta_path[active_set, k-1] != 0), ][which(sub_mismatch == 1)]
      delta_rho = min(delta_rho_sub_tmp)
      delta_rho_vec = c(delta_rho_vec, delta_rho)
      
      # predictor shrink to zero
      sub_viol = which(active_set)[which(beta_path[which(active_set), k-1] / delta_b == delta_rho)]
      
    } else {
      delta_rho_vec = c(delta_rho_vec, -1) 
    }
    
    # 2. dual feasibility & KKT - stationarity condition
      # for NOT active set
    if(sum(active_set) != p) { # if all predictors are activated, do not run opt problem(all predictors are actived)
      delta = Variable(1)
      constraints = list(-t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta * delta_b)) + 
                           t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta * delta_l) 
                         <= (rho_path[k-1] - delta) * matrix(rep(1, sum(!active_set)), ncol = 1),
                         
                         -t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta * delta_b)) +
                           t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta * delta_l) 
                         >= -(rho_path[k-1] - delta) * matrix(rep(1, sum(!active_set)), ncol = 1),
                         
                         delta >= 0)
      objective = Maximize(delta)
      problem = Problem(objective, constraints)
      result = solve(problem)
      
      # find predictor which cannot satisfy KKT condition among NOT active set
        # -> this must be actived
      kkt_viol_to_active = c()
      
      # NA for delta value : there exist some predictor which never satisfy KKT condition
      if(is.na(result$getValue(delta))) { 
        delta_rho_vec = c(delta_rho_vec, -1)
        
        # predictors which cannot satisfy KKT condition
        kkt_viol_to_active1 = which(-t(X[, !active_set]) %*% (y - X[, active_set] %*% beta_path[active_set, k-1]) + 
                                      t(Aeq[, !active_set]) %*% lambda_patheq[, k-1] >
                                      rho_path[k-1] * matrix(rep(1, sum(!active_set)), ncol = 1))
        kkt_viol_to_active2 = which(-t(X[, !active_set]) %*% (y - X[, active_set] %*% beta_path[active_set, k-1]) + 
                                      t(Aeq[, !active_set]) %*% lambda_patheq[, k-1] <
                                      - rho_path[k-1] * matrix(rep(1, sum(!active_set)), ncol = 1))
        kkt_viol_to_active = sort(union(kkt_viol_to_active1, kkt_viol_to_active2)) # next added predictor
        
      } else {
        delta_rho_vec = c(delta_rho_vec, result$getValue(delta))
        
        # predictor on boundary(for update active set) -> update for active_set
        delta_rho_kkt_tmp = result$getValue(delta)
        kkt_viol1 = which(abs(-t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta_rho_kkt_tmp * delta_b)) + 
                                t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta_rho_kkt_tmp * delta_l) - 
                                (rho_path[k-1] - delta_rho_kkt_tmp) * matrix(rep(1, sum(!active_set)), ncol = 1)) <= 1e-6)
        kkt_viol2 = which(abs(-t(X[, !active_set]) %*% (y - X[, active_set] %*% (beta_path[active_set, k-1] - delta_rho_kkt_tmp * delta_b)) + 
                                t(Aeq[, !active_set]) %*% (lambda_patheq[, k-1] - delta_rho_kkt_tmp * delta_l) + 
                                (rho_path[k-1] - delta_rho_kkt_tmp) * matrix(rep(1, sum(!active_set)), ncol = 1)) <= 1e-6)
        kkt_viol = sort(union(kkt_viol1, kkt_viol2)) # next added predictor
      }
      
      # all predictors are activated
    } else { 
      kkt_viol_to_active = c()
      delta_rho_vec = c(delta_rho_vec, -1)
    }
    
    ### terminate checking (rho becomes 0)
    if(any(delta_rho_vec > 0)) {
      if(min(delta_rho_vec[delta_rho_vec > 0]) > rho_path[k-1]) {
        delta_rho = rho_path[k-1]
        rho_path[k] = rho_path[k-1] - delta_rho
        beta_path[active_set, k] = beta_path[active_set, k-1, drop=F] - delta_rho * delta_b
        objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
        lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_l
        break
      }
    }
    
    ### choose delta_rho
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
      
      # 4. there is no candidates for delta_rho 
    } else if(delta_rho_vec[1] < 0 & delta_rho_vec[2] < 0){
      
      # if we need to change active set (NOT terminate algorithm)
      if(length(kkt_viol_to_active) > 0) { 
        rho_path[k] = rho_path[k-1]
        beta_path[active_set, k] = beta_path[active_set, k-1, drop=F]
        objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
        lambda_patheq[, k] = lambda_patheq[, k-1]
        
        active_set[which(!active_set)[kkt_viol_to_active]] = T
        num_active = sum(active_set)
        next
        
      } else {
        ### terminate algorithm (no event is upcoming)
        delta_rho = rho_path[k-1] # make rho to zero
        rho_path[k] = rho_path[k-1] - delta_rho
        beta_path[active_set, k] = beta_path[active_set, k-1, drop=F] - delta_rho * delta_b
        objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
        lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_l
        break
      }
    }
    
    ### updates
    rho_path[k] = rho_path[k-1] - delta_rho
    beta_path[active_set, k] = beta_path[active_set, k-1, drop=F] - delta_rho * delta_b
    objval_path[k] = sum((y - X %*% beta_path[ ,k])^2) + rho_path[k] * sum(abs(beta_path[k]))
    lambda_patheq[, k] = lambda_patheq[, k-1] - delta_rho * delta_l
    
      # case 1) subgradient rule is violated
    if(min_idx == 1) {
      active_set[sub_viol] = F
      
      # case 2) KKT condition is violated
    } else if(min_idx == 2) {
      active_set[which(!active_set)[kkt_viol]] = T
    } 
    
    if(length(kkt_viol_to_active) > 0) {
      active_set[which(!active_set)[kkt_viol_to_active]] = T
    }
    num_active = sum(active_set)
  }
  
  if(show_plotting == T) {
    # plot
    # 1. by sum(abs(beta))
    par(mfrow=c(1,1))
    plot(apply(abs(beta_path[, 1:k]), 2, sum), rep(0, k), type = "l", col = 1, lwd = 1, 
         xlim = c(0, sum(abs(beta_path[, k]))), ylim = c(min(beta_path[, 1:k]), max(beta_path[, 1:k])),
         xlab = "sum(abs(beta))", ylab = "beta coef", main = "beta coef path(by sum(abs(beta)))")
    for(i in 1:p) lines(apply(abs(beta_path[, 1:k]), 2, sum), beta_path[i, 1:k], type = "l", col = i+1, lwd = 2)
    
    # 2. by step
    par(mfrow=c(1,1))
    plot(seq(1, k), rep(0, k), type = "l", col = 1, lwd = 1, 
         xlim = c(1, k), ylim = c(min(beta_path[, 1:k]), max(beta_path[, 1:k])),
         xlab = "steps", ylab = "beta coef", main = "beta coef path(by step)")
    for(i in 1:p) lines(seq(1, k), beta_path[i, 1:k], type = "l", col = i+1, lwd = 2)
  }
  
  # equality constraint check
  flag = (sum(abs(Aeq %*% beta_path[, k])) <= 1e-6)
  
  # result
  return(list(rho_path = rho_path[1:k],
              beta_path = beta_path[, 1:k],
              lambda_patheq = lambda_patheq[, 1:k],
              objval_path = objval_path[1:k],
              step_num = k,
              constraint_check = flag))
}