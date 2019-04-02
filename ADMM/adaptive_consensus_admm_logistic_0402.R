### Adaptive Consensus ADMM for Distributed Optimization (4/2)###
if(!require(data.table)) install.packages('data.table')
require(data.table)

# 1. parameter setting
k = 10 # number of machines
n = 5000
p = 100
eta = 0.05
### penalty parameters for each slaves
penalty_vec = rep(0.001, k)
T_f = 2
e_cor = 0.2
e_tol = 0.1
C_cg = 10e+10
epochs = 100
global_objs = list(k = k, n = n, p = p, eta = eta, penalty_vec = penalty_vec)
### true beta
true_beta = c(rep(1, 5), rep(0, p - 5))
set.seed(520)
# global_beta = matrix(rnorm(p, mean = 0, sd = 0.1), p, 1)
global_beta = matrix(rep(0, p), p, 1)
### local beta, dual
local_beta = list()
local_dual = list()
local_dual_hat = list()
for(i in 1:k)
{
  set.seed(i * 520)
  # local_beta[[i]] = rnorm(p, mean = 0, sd = 0.1)
  local_beta[[i]] = matrix(rep(0, p), p, 1)
  set.seed(i * 960)
  # local_dual[[i]] = rnorm(p, mean = 0, sd = 0.1)
  local_dual[[i]] = matrix(rep(0, p), p, 1)
  local_dual_hat[[i]] = rep(0, p)
}

# 2. data generation
for(i in 1:k)
{
  set.seed(i)
  x = matrix(rnorm((n / k) * p, i, i * 0.1), (n / k), p)
  y = rbinom((n / k), size = 1, prob = 1 / (1 + exp(-x %*% true_beta)))
  fwrite(cbind(x, y), file = paste('C:/Users/dpelt/r_default/admm_adaptive_data/data', i, '.csv', sep=""))
}

# 3. stopping crietria
pri_crieria_fun = function(global_beta, local_beta)
{
  local_beta_norm = 0
  for(i in 1:k) local_beta_norm = local_beta_norm + norm(as.matrix(local_beta[[i]]), "F")
  return(e_tol * max(local_beta_norm, k * norm(global_beta, "F")))
}
dual_criteria_fun = function(local_dual)
{
  local_dual_norm = 0
  for(i in 1:k) local_dual_norm = local_dual_norm + norm(as.matrix(local_dual[[i]]), "F")
  return(e_tol * local_dual_norm)
}

# 4. slave function for grad computing
slave_grad_func = function(index)
{
  data = as.matrix(fread(file = paste('C:/Users/dpelt/r_default/admm_adaptive_data/data', index, '.csv', sep=""), header = T))
  x = data[, 1:global_objs$p]
  y = as.matrix(data[, (global_objs$p + 1)])
  b = local_beta[[index]]
  dual = local_dual[[index]]
  # while(T)
  # {
  #   # computing gradient
  #   grad = -matrix((t(y) %*% x) / (global_objs$n / k), global_objs$p, 1) + (t(x) %*% (1 / (1 + exp(-x %*% b)))) / (global_objs$n / k) +
  #          penalty_vec[index] * matrix(global_beta - b + dual / penalty_vec[index], global_objs$p, 1)
  #   # update local beta
  #   b = b - global_objs$eta * grad
  #   # checking stopping criteria
  #   if(norm(grad, "F") < 1e-6) break
  # }
  cat(index, "\n")
  for(j in 1:5000)
  {
    # computing gradient
    grad1 = -matrix((t(y) %*% x) / (global_objs$n / k), global_objs$p, 1) + (t(x) %*% (1 / (1 + exp(-x %*% b)))) / (global_objs$n / k)
    grad2 = global_objs$penalty_vec[index] * matrix(global_beta - b + (dual / global_objs$penalty_vec[index]), global_objs$p, 1)
    # update local beta
    b = b - (global_objs$eta / (n / k)) * (grad1 + grad2)
    # checking stopping criteria
    if(j %% 500 == 0) cat(norm(grad1, "F"), norm(grad2, "F"), "\n")
    # if(norm(grad, "F") < 0.1) break
  }
  return(b)
}

# 5. slave function for penalty parameter update
slave_penalty_func = function(i)
{
  updated_local_dual_hat = local_dual[[i]] + global_objs$penalty_vec[i] * (global_beta - updated_local_beta[[i]])
  grad_local_beta = updated_local_beta[[i]] - local_beta[[i]]
  grad_local_dual = updated_local_dual[[i]] - local_dual[[i]]
  grad_local_dual_hat = updated_local_dual_hat - local_dual_hat[[i]]
  grad_global_beta = global_beta - updated_global_beta
  ###
  alpha_SD = (t(grad_local_dual_hat) %*% grad_local_dual_hat) / (t(grad_local_beta) %*% grad_local_dual_hat)
  alpha_MG = (t(grad_local_beta) %*% grad_local_dual_hat) / (t(grad_local_beta) %*% grad_local_beta)
  beta_SD = (t(grad_local_dual) %*% (grad_local_dual)) / (t(grad_global_beta) %*% grad_local_dual)
  beta_MG = (t(grad_global_beta) %*% grad_local_dual_hat) / (t(grad_global_beta) %*% grad_global_beta)
  ###
  if(2 * alpha_MG > alpha_SD) {
    alpha = alpha_MG
  } else {
    alpha = (alpha_SD - alpha_MG / 2)
  }
  if(2 * beta_MG > beta_SD) {
    beta = beta_MG
  } else {
    beta = (beta_SD - beta_MG / 2)
  }
  ###
  alpha_corr = (t(grad_local_beta) %*% grad_local_dual_hat) / (norm(grad_local_beta, "F") * norm(grad_local_dual_hat, "F"))
  beta_corr = (t(grad_global_beta) %*% grad_local_dual) / (norm(grad_global_beta, "F") * norm(grad_local_dual, "F"))
  if(alpha_corr > e_cor && beta_corr > e_cor) {
    penalty_hat = 1 / sqrt(alpha * beta)
  } else if(alpha_corr > e_cor && beta_corr <= e_cor) {
    penalty_hat = 1 / alpha
  } else if(alpha_corr <= e_cor && beta_corr > e_cor) {
    penalty_hat = 1 / beta
  } else penalty_hat = global_objs$penalty_vec[i]
  penalty = max(min(penalty_hat, (1 + C_cg / t^2) * global_objs$penalty_vec[i]), global_objs$penalty_vec[i] / (1 + C_cg / t^2))
  ###
  result = list(updated_local_dual_hat = updated_local_dual_hat, penalty = penalty)
  return(result)
}

# 6. Adaptive consensus ADMM(this would be upgraded with mpi)
for(t in 1:epochs)
{
  ### update local beta
  updated_local_beta = list()
  for(i in 1:k) updated_local_beta[[i]] = slave_grad_func(i)
  ### update global beta
  sum = rep(0, p)
  for(i in 1:k) sum = sum + updated_local_beta[[i]] - local_dual[[i]] / global_objs$penalty_vec
  updated_global_beta = sum / sum(global_objs$penalty_vec)
  fwrite(as.matrix(updated_global_beta), file = paste('C:/Users/dpelt/r_default/admm_adaptive_beta/beta', t, '.csv', sep=""))
  ### update local dual
  updated_local_dual = list()
  for(i in 1:k) updated_local_dual[[i]] = local_dual[[i]] + global_objs$penalty_vec[i] * (updated_global_beta - updated_local_beta[[i]])
  ### update penalty parameter
  if(t %% T_f == 0)
  {
    for(i in 1:k) 
    {
      result = slave_penalty_func(i)
      global_objs$penalty_vec[i] = result$penalty
      local_dual_hat[[i]] = result$updated_local_dual_hat
    }
  } else {
    for(i in 1:k) local_dual_hat[[i]] = local_dual[[i]] + global_objs$penalty_vec[i] * (global_beta - updated_local_beta[[i]])
  }
  ###
  global_beta = updated_global_beta
  local_beta = updated_local_beta
  local_dual = updated_local_dual
}

# 7. integrating global beta
total_beta = matrix(numeric(0), p, 0)
for(i in 1:t)
{
  total_beta = cbind(total_beta, as.matrix(fread(file = paste('C:/Users/dpelt/r_default/admm_adaptive_beta/beta', i, '.csv', sep=""))))
}

# 8. loss computing
loss = rep(0, t)
for(i in 1:k)
{
  data = as.matrix(fread(file = paste('C:/Users/dpelt/r_default/admm_adaptive_data/data', i, '.csv', sep=""), header = T))
  x = data[ ,1:p]
  y = data[ ,(p + 1)]
  loss = loss - t(y) %*% x %*% total_beta + colSums(log(1 + exp(x %*% total_beta)))
}
fwrite(as.matrix(loss), file = paste('C:/Users/dpelt/r_default/adaptive_consensus_loss.csv', sep=""))

# 8. save loss~running_time plot
# png(file = "'C:/Users/dpelt/r_default/adaptive_consensus_loss_plot.png")
plot(seq(1, t, 1), loss)
# dev.off()