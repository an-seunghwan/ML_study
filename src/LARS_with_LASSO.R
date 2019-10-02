rm(list=ls())
gc()
set.seed(520)
if(!require("lars")) install.packages('lars')
library('lars')

### LARS with LASSO algorithm ###
# use longley data #

# setting
raw_X = data.matrix(longley[, 1:6])
  # check linearly independent of data matrix
if(qr(raw_X)$rank != ncol(raw_X)) {
  cat("Data matrix is not Linearly Independent!\n")
  return()
}
X_scaled = cbind(rep(1, 16), scale(raw_X)) # include intercept
n = dim(X_scaled)[1]
p = dim(X_scaled)[2]
y = data.matrix(longley[, 7])
y_scaled = y - mean(y)

# LARS function
lars_lasso_function = function(X_scaled, y_scaled, plotting = F) {
  ### functions
  get_alpha_plus = function(r, x, new_x, X, d) {
    alpha_plus = (t(r) %*% x - t(r) %*% new_x) / (t(r) %*% x - t(X %*% d) %*% new_x)
    return(alpha_plus)
  }
  
  get_alpha_minus = function(r, x, new_x, X, d) {
    alpha_minus = (t(r) %*% x + t(r) %*% new_x) / (t(r) %*% x + t(X %*% d) %*% new_x)
    return(alpha_minus)
  }
  
  get_alpha_zero = function(beta, d) {
    alpha = - beta / d
    return(alpha)
  }
  
  index_set = c() # 현재 step까지 선택된 predictor들의 인덱스를 저장
  beta_path = matrix(numeric(0), nrow = p) # beta의 결과 저장
  
  ### init
  r = y_scaled
  beta = rep(0, p)
  beta_path = cbind(beta_path, beta)
  index = which.max(abs(t(X_scaled) %*% r)) 
  index_set = c(index_set, index)
  E = matrix(numeric(0), nrow = p)
  flag = 1 # alpha에 대한 case 확인
  
  ### 모든 predictor가 선택될 때까지 진행
  while(length(unique(index_set)) < p) {
    if(flag == 1) {
      e = rep(0, p)
      e[index] = 1
      E = cbind(E, e)
    } else if(flag == 2) {
      E = E[, -remove_id]
    }
    
    # direction
    d = E %*% solve(t(X_scaled %*% E) %*% X_scaled %*% E) %*% t(X_scaled %*% E) %*% r
    
    # alpha
    # first case
    NOT_index_set = seq(1, p, 1)[!(seq(1, p, 1) %in% index_set)]
    alpha_vec1 = c()
    for(i in 1:length(NOT_index_set)) {
      alpha_plus = get_alpha_plus(r, X_scaled[, index], X_scaled[, NOT_index_set[i]], X_scaled, d)
      alpha_minus = get_alpha_minus(r, X_scaled[, index], X_scaled[, NOT_index_set[i]], X_scaled, d)
      alpha_vec1 = c(alpha_vec1, alpha_plus, alpha_minus)
    }
    
    # second case
    alpha_vec2 = c()
    for(i in index_set[index_set != index]) {
      alpha = get_alpha_zero(beta[i], d[i])
      alpha_vec2 = c(alpha_vec2, alpha)
    }
    
    # choice
    alpha1 = min(alpha_vec1[0 < alpha_vec1 & alpha_vec1 < 1])
    if(sum(0 < alpha_vec2 & alpha_vec2 < 1) > 0) {
      alpha2 = min(alpha_vec2[0 < alpha_vec2 & alpha_vec2 < 1])
    } else alpha2 = 1
    
    # 2 cases 
    if(alpha1 < alpha2) {
      index = NOT_index_set[ifelse(which(alpha_vec1 == alpha1) %% 2 == 0, which(alpha_vec1 == alpha1) / 2, 
                                   (which(alpha_vec1 == alpha1) + 1) / 2)]
      index_set = c(index_set, index)
      alpha = alpha1
      flag = 1
    } else {
      index = index_set[index_set != index][which(alpha_vec2 == alpha2)]
      remove_id = which(index_set == index)
      index_set = index_set[index_set != index]
      alpha = alpha2
      flag = 2
    }
    
    # update
    beta = beta + alpha * d
    r = y_scaled - X_scaled %*% beta
    beta_path = cbind(beta_path, beta)
  }
  
  ### for the last step
  if(flag == 1) {
    e = rep(0, p)
    e[index] = 1
    E = cbind(E, e)
  } else if(flag == 2) {
    E = E[, -index]
  }
  
  # direction
  d = E %*% solve(t(X_scaled %*% E) %*% X_scaled %*% E) %*% t(X_scaled %*% E) %*% r
  
  # alpha
  alpha = 1
  
  # update
  beta = beta + alpha * d
  r = y_scaled - X_scaled %*% beta
  beta_path = cbind(beta_path, beta)
  
  # plotting
  if(plotting) {
    par(mfrow = c(1,1))
    beta_max = max(beta_path)
    beta_min = min(beta_path)
    sum_abs_beta = apply(abs(beta_path), 2, sum)
    plot(sum_abs_beta, beta_path[1, ], ylim = c(beta_min, beta_max), type = 'l', col = 1, lwd = 2,
         xlab = "sum(abs(beta))", ylab = "beta", main = 'BETA coef PATH(LASSO)')
    for(i in 2:p) points(sum_abs_beta, beta_path[i, ], type = 'l', col = i, lwd = 2)
  }
  
  return(beta_path)
}

# result
beta_path = lars_lasso_function(X_scaled, y_scaled, plotting = T)
step = dim(beta_path)[2]

# check
my_beta = beta_path[, 12]
lars_beta = coef(lars(X_scaled, y_scaled, type="lar"))[p, ]
cat("step num: ", step, "\n")
cat("my beta: ", my_beta, "\n")
cat("lars beta: ", lars_beta, "\n")
