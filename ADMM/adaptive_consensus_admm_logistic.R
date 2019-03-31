### Adaptive Consensus ADMM for Distributed Optimization ###
if(!require(data.table)) install.packages('data.table')
require(data.table)

# 1. parameter setting
k = 10 # number of machines
n = 5000
p = 100
eta = 0.1
## penalty parameters
penalty_vec = rep(1, k)
T_f = 2
e_cor = 0.2
e_tol = 0.1
C_cg = 10e+10
epochs = 100
global_objs = list(k = k, n = n, p = p, eta = eta, penalty_vec = penalty_vec)
## true beta
idx = sample(1 : p, size = 5, replace = F)
true_beta = rep(0, p)
true_beta[idx] = 1
## global beta
set.seed(520)
global_beta = matrix(rnorm(p), p, 1)
## local beta, dual
local_beta = list()
local_dual = list()
for(i in 1:k)
{
  set.seed(i * 500)
  local_beta[[i]] = rnorm(p)
  set.seed(i * 960)
  local_dual[[i]] = rnorm(p)
}
old_dual_hat = list()

# 2. data generation
for(i in 1:k)
{
  set.seed(i)
  x = matrix(rnorm((n / k) * p, (i %% 3), i), (n / k), p)
  y = rbinom((n / k), size = 1, prob = 1 / (1 + exp(-x %*% true_beta)))
  fwrite(cbind(x, y), file = paste('C:/Users/dpelt/r_default/admm_adaptive_data/data', i, '.csv', sep=""))
}

# 3. stopping crietria
pri_crieria_fun = function(global_beta, local_beta)
{
  local_beta_norm = 0
  for(i in 1:k)
  {
    local_beta_norm = local_beta_norm + norm(as.matrix(local_beta[[i]]), "F")
  }
  return(e_tol * max(local_beta_norm, k * norm(global_beta, "F")))
}
dual_criteria_fun = function(local_dual)
{
  local_dual_norm = 0
  for(i in 1:k)
  {
    local_dual_norm = local_dual_norm + norm(as.matrix(local_dual[[i]]), "F")
  }
  return(e_tol * local_dual_norm)
}

# 4. slave function
slavefunction = function(index)
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
  #     penalty_vec[index] * matrix(global_beta - b + dual / penalty_vec[index], global_objs$p, 1)
  #   # update local beta
  #   b = b - global_objs$eta * grad
  #   # checking stopping criteria
  #   if(norm(grad, "F") < 1e-6) break
  # }
  for(j in 1:1000)
  {
    # computing gradient
      grad = -matrix((t(y) %*% x) / (global_objs$n / k), global_objs$p, 1) + (t(x) %*% (1 / (1 + exp(-x %*% b)))) / (global_objs$n / k) +
        global_objs$penalty_vec[index] * matrix(-global_beta + b + dual / global_objs$penalty_vec[index], global_objs$p, 1)
    # update local beta
    b = b - global_objs$eta * grad
    # checking stopping criteria
    # cat(norm(grad, "F"), "\n")
    # if(norm(grad, "F") < 0.1) break
  }
  return(b)
}

# 5. Adaptive consensus ADMM(this would be upgraded with mpi)
for(t in 1:epochs)
{
  # update local beta
  updated_local_beta = list()
  for(i in 1:k)
  {
    updated_local_beta[[i]] = slavefunction(i)
  }
  # update global beta
  sum = rep(0, p)
  for(i in 1:k) sum = sum + updated_local_beta[[i]] + local_dual[[i]] / global_objs$penalty_vec[i]
  old_beta = global_beta
  global_beta = sum / k
  fwrite(as.matrix(global_beta), file = paste('C:/Users/dpelt/r_default/admm_adaptive_beta/beta', t, '.csv', sep=""))
  # update local dual
  updated_local_dual = list()
  for(i in 1:k) updated_local_dual[[i]] = local_dual[[i]] + global_objs$penalty_vec[i] * (-global_beta + updated_local_beta[[i]])
  if(t %% T_f == 0)
  {
    dual_hat = list()
    alpha_corr = rep(0, k)
    beta_corr = rep(0, k)
    penalty_hat = rep(0, k)
    for(i in 1:k)
    {
      dual_hat[[i]] = local_dual[[i]] + global_objs$penalty_vec[i] * (old_beta - local_beta[[i]])
      fit1 = lm(updated_local_beta[[i]] ~ dual_hat[[i]])
      alpha = fit1$coefficients[[2]]
      fit2 = lm(global_beta ~ local_dual[[i]])
      beta = fit2$coefficients[[2]]
      alpha_corr[i] = t(updated_local_beta[[i]] - local_beta[[i]]) %*% (dual_hat[[i]] - old_dual_hat[[i]]) / norm(updated_local_beta[[i]] - local_beta[[i]], "F") * norm(dual_hat[[i]] - old_dual_hat[[i]], "F")
      beta_corr[i] = t(-global_beta + old_beta) %*% (updated_local_dual[[i]] - local_dual[[i]]) / norm(-global_beta + old_beta, "F") * norm(updated_local_dual[[i]] - local_dual[[i]], "F")
      if(alpha_corr[i] > e_cor && beta_corr[i] > e_cor) {
        penalty_hat[i] = 1 / sqrt(alpha * beta)
      } else if(alpha_corr[i] > e_cor && beta_corr[i] <= e_cor) {
        penalty_hat[i] = 1 / alpha
      } else if(alpha_corr[i] <= e_cor && beta_corr[i] > e_cor) {
        penalty_hat[i] = 1 / beta
      } else global_objs$penalty_hat[i] = global_objs$penalty_vec[i]
      global_objs$penalty_vec[i] = max(min(penalty_hat[i], (1 + C_cg / t^2) * global_objs$penalty_vec[i]), global_objs$penalty_vec[i] / (1 + C_cg / t^2))
    }
  } else
  {
    for(i in 1:k)
    {
      old_dual_hat[[i]] = local_dual[[i]] + global_objs$penalty_vec[i] * (old_beta - local_beta[[i]])
    }
  }
  for(i in 1:k)
  {
    local_beta[[i]] = updated_local_beta[[i]]
    local_dual[[i]] = updated_local_dual[[i]]
  }
}

# 6. integrating global beta
total_beta = matrix(numeric(0), p, 0)
for(i in 1:t)
{
  total_beta = cbind(total_beta, as.matrix(fread(file = paste('C:/Users/dpelt/r_default/admm_adaptive_beta/beta', i, '.csv', sep=""))))
}

# 7. loss computing
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
plot(seq(1,t,1), loss)
# dev.off()