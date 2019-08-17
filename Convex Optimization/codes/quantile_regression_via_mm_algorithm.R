rm(list=ls())
gc()

# quantile regression example

# 1. paramters
set.seed(520)
n = 100
p = 10
q = 1/2
threshold = 1e-6
true_b = matrix(runif(10, min = 1, max = 3), nrow = p)
X = matrix(rnorm(n * p), nrow = n, byrow = T)
e = matrix(rnorm(n, 0, 0.3), nrow = n)
y = X %*% true_b + e

# 2. simulation
beta = matrix(rep(1, p), nrow = p)
for(iter in 1:100) {
  r = c(y - X %*% beta)
  W = diag(1 / (4 * abs(r)))
  term1 = solve(t(X) %*% W %*% X)
  term2 = t(X) %*% W %*% y
  term3 = colSums(X)
  
  new_beta = term1 %*% term2 + (q / 2 - 1 / 4) * term1 %*% term3
  if(max(abs(beta - new_beta)) < threshold) break 
  beta = new_beta # update
}
print(new_beta)
print(beta)
