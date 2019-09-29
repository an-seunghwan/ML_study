#
init_OPL_lasso = function(P, Aeq, beq)
{
  ################### initialization #####################
  p = dim(P)[1]
  m = dim(Aeq)[1]
  
  zero_m = matrix(rep(0, m), nrow = m)
  zero_p = matrix(rep(0, p), nrow = p)
  one_p = matrix(rep(1, p), nrow = p)
  one_m = matrix(rep(1, m), nrow = m)
  I_pp = diag(1, p)
  zero_pm = matrix(rep(0, p*m), nrow = p)
  zero_pp = matrix(rep(0, p*p), nrow = p)
  zero_mm = matrix(rep(0, m*m), nrow = m)
  
  ### get initial x (x_init)
  ##### solve
  # min   sum(abs(beta))
  # s.t.  Ax = b
  #       x >= 0
  #####
  x = Variable(p)
  constraints = list(Aeq %*% x == beq,
                     x >= 0)
  objective = Minimize(sum(abs(x)))
  problem = Problem(objective, constraints)
  result = solve(problem)
  x_init = result$getValue(x)
  x_init[x_init < 1e-4] = 0.0 # for numerical stability
  x_init
  
  # ### conditions check
  # sum(abs(x_init))
  # abs(Aeq %*% x_init - beq) < 1e-6
  # x_init >= 0
  
  ### get initial lambda, mu, nu, active_set
  # Design matrix
  D1 = cbind(1, t(zero_m), t(zero_p))
  D2 = cbind(zero_p, t(Aeq), I_pp)
  D3 = cbind(one_p, zero_pm, zero_pp)
  D4 = cbind(zero_p, zero_pm, I_pp)
  
  # get initial lambda, mu, nu corresponding to x_init
  target = Variable(1 + m + p)
  constraints = list(2 * P %*% x_init + D2 %*% target <= D3 %*% target, # stationarity
                     2 * P %*% x_init + D2 %*% target >= - D3 %*% target, # stationarity
                     D1 %*% target >= 0, # lambda must be positive
                     D4 %*% target >= zero_p, # dual feasibility
                     t(D4 %*% target) %*% x_init == 0 # complementary slackness
                     )
  objective = Minimize(D1 %*% target)
  problem = Problem(objective, constraints)
  result = solve(problem)
  # result$getValue(target)
  max_lambda = result$getValue(target)[1]
  max_mu = result$getValue(target)[2:(1+m), , drop=F]
  max_nu = result$getValue(target)[(m+2):nrow(result$getValue(target)), , drop=F]
  max_nu[max_nu < 1e-4] = 0.0 # for numerical stability 
  # x_init; max_lambda; max_mu; max_nu
  
  # active_set
  active_set = rep(F, p)
  active_set[which(x_init != 0)] = T
  # active_set
  # sum(active_set)
  
  return(list(max_lambda = max_lambda,
              x_init = x_init,
              max_mu = max_mu,
              max_nu = max_nu,
              active_set = active_set))
}