rm(list=ls())
gc()
set.seed(520)

# coordinatewise alogrithm for lasso regression

# 1. parameters
n = 1000
p = 10
true_b = c(1, 2, 1, 3, 2, 1, 3, 2, 1, 3)
X = matrix(rnorm(n * p), nrow = n, byrow = F)
e = matrix(rnorm(n), nrow = n)
y = X %*% true_b + e
lambda = 1
max_iter = 1000
threshold = 1e-6

#
sign = function(x) {
  if(x > 0) return(1)
  else if(x < 0) return(-1)
  else return(0)
}

# 2. coordinatewise alogrithm
b = rnorm(p)
while(T) {
  new_b = b
  for(index in 1:p) {
    new_b[index] = 0
    r_term = y - X %*% new_b
    term1 = t(r_term) %*% X[, index]
    term2 = sum(X[, index] ^ 2)
    if(abs(-2 * term1) > lambda) new_b_index = term1 / term2 + sign(-2 * term1) * lambda / (2 * term2)
    else new_b_index = 0
    new_b[index] = new_b_index
  }
  if(max(abs(new_b - b)) < threshold) break
  b = new_b
}

# 3. result
estb = b
cat("estimate: ", estb, "\n")
cat("true: ", true_b, "\n")