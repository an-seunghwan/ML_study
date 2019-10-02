rm(list=ls())
gc()
set.seed(520)
if(!require("lars")) install.packages('lars')
library('lars')

### LARS algorithm ###
# 1. where n > p
# 2. NOT lasso version, only regression
# 3. for scaled, linearly independent data

# setting
n = 1000
p = 10
X = matrix(rnorm(n*p, mean = 1, sd = 3), nrow = n)
X_scaled = scale(X)
y = rnorm(n, mean = 1, sd = 3)
y_scaled = y - mean(y)

# LARS function
lars_function = function(X_scaled, y_scaled, plotting = F) {
  # check linear independence
  if(qr(X)$rank != ncol(X)) {
    cat("Data matrix is not Linearly Independent!\n")
    return()
  }
  
  get_alpha_plus = function(r, x, new_x, X, d) {
    alpha_plus = (t(r) %*% x - t(r) %*% new_x) / (t(r) %*% x - t(X %*% d) %*% new_x)
    return(alpha_plus)
  }
  
  get_alpha_minus = function(r, x, new_x, X, d) {
    alpha_minus = (t(r) %*% x + t(r) %*% new_x) / (t(r) %*% x + t(X %*% d) %*% new_x)
    return(alpha_minus)
  }
  index_set = c() # 현재 step까지 선택된 predictor들의 인덱스를 저장
  beta_path = matrix(numeric(0), nrow = p) # beta의 결과 저장
  
  # init
  r = y_scaled
  beta = rep(0, p)
  beta_path = cbind(beta_path, beta)
  index = which.max(abs(t(X_scaled) %*% r)) # index means lastest chosen index
  index_set = c(index_set, index)
  E = matrix(numeric(0), nrow = p)
  
  # LARS
  for(step in 1:p) { 
    if(step != p){
      e = rep(0, p)
      e[index] = 1
      E = cbind(E, e)
      
      # direction
      d = E %*% solve(t(X_scaled %*% E) %*% X_scaled %*% E) %*% t(X_scaled %*% E) %*% r
      
      # alpha
      NOT_index_set = seq(1, p, 1)[!(seq(1, p, 1) %in% index_set)]
      alpha_vec = c()
      for(i in 1:length(NOT_index_set)) {
        alpha_plus = get_alpha_plus(r, X_scaled[, index], X_scaled[, NOT_index_set[i]], X_scaled, d)
        alpha_minus = get_alpha_minus(r, X_scaled[, index], X_scaled[, NOT_index_set[i]], X_scaled, d)
        alpha_vec = c(alpha_vec, alpha_plus, alpha_minus)
      }
      alpha = min(alpha_vec[0 < alpha_vec & alpha_vec < 1])
      index = NOT_index_set[ifelse(which(alpha_vec == alpha) %% 2 == 0, which(alpha_vec == alpha) / 2, 
                                   (which(alpha_vec == alpha) + 1) / 2)]
      
      # update
      beta = beta + alpha * d
      r = y_scaled - X_scaled %*% beta
      index_set = c(index_set, index)
      beta_path = cbind(beta_path, beta)
    }
    else { # for last step(p)
      e = rep(0, p)
      e[index] = 1
      E = cbind(E, e)
      
      # direction
      d = E %*% solve(t(X_scaled %*% E) %*% X_scaled %*% E) %*% t(X_scaled %*% E) %*% r
      
      # alpha
      alpha = 1
      
      # update
      beta = beta + alpha * d
      r = y_scaled - X_scaled %*% beta
      beta_path = cbind(beta_path, beta)
    }
  }
  
  # plotting
  if(plotting) {
    par(mfrow = c(1,1))
    beta_max = max(beta_path)
    beta_min = min(beta_path)
    plot(beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
         xlab = "step + 1", ylab = "beta", main = 'beta coefficients path')
    for(i in 2:p) points(beta_path[i, ], type = 'l', col = i, lwd = 2)
  }
  
  return(beta_path)
}

# result
beta_path = lars_function(X_scaled, y_scaled, plotting = T)

# check
my_beta = beta_path[, p + 1]
lm_beta = lm(y_scaled ~ X_scaled - 1)$coefficients
lars_beta = coef(lars(X_scaled, y_scaled, type="lar"))[p + 1, ]
cat("my beta: ", my_beta, "\n")
cat("lm beta: ", lm_beta, "\n")
cat("lars beta: ", lars_beta, "\n")