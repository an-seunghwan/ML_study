rm(list=ls())
gc()

# logistic regression via MM algorithm

# 1. parameters
set.seed(520)
n = 100
p = 3
x1 = rnorm(n, sd = 0.8)
x2 = rnorm(n, sd = 0.6)
true_b = c(1, 2, 3)
z = true_b[1] + true_b[2] * x1 + true_b[3] * x2
prob = exp(z) / (1 + exp(z))
y = rbinom(n, 1, prob)
X = cbind(rep(1, n), x1, x2)
df = data.frame(y = y, x1 = x1, x2 = x2)
max_iter = 100000
threshold = 1e-6

# 2. MM alogrithm
b = matrix(rep(1, 3), nrow = p)
term = solve(t(X) %*% X)
start1 = Sys.time()
for(iter in 1:max_iter) {
  theta = exp(X %*% b) / (1 + exp(X %*% b))
  new_b = b - 4 * term %*% t(X) %*% (theta - y)
  if(max(abs(new_b - b)) < threshold) break
  b = new_b
}
end1 = Sys.time()
estb_MM = c(b)

# 3. NR
b = matrix(rep(1, 3), nrow = p)
start2 = Sys.time()
for(iter in 1:max_iter){
  theta = exp(X %*% b)/(1 + exp(X %*% b))
  grad = -t(X) %*% (y - theta)
  hess = drop(t(theta)%*%(1 - theta)) * t(X) %*% X
  new_b = b - solve(hess) %*% grad
  if (max(abs(new_b - b)) < threshold) break
  b = new_b  
}
end2 = Sys.time()
estb_NR = c(b)

# 4. glm
estb_glm = glm(y ~ x1 + x2, data = df, family = "binomial")

# 5. l2 package
b = matrix(rep(1, 3), nrow = p)
start3 = Sys.time()
for(iter in 1:max_iter) {
  theta = exp(X %*% b)/(1 + exp(X %*% b))
  tilde_X = X / 2
  tilde_y = 2 * (y - theta) - X %*% b / 2
  new_b = lm(tilde_y ~ tilde_X - 1)$coefficient # l2 package...?
  if(max(abs(new_b - b)) < threshold) break
  b = new_b
}
end3 = Sys.time()
estb_l2 = c(b)
# Warning message:
  # In max(abs(new_b - b)) : no non-missing arguments to max; returning -Inf

# results
estb_glm$coefficients
estb_MM ; end1 - start1
estb_NR ; end2 - start2
estb_l2 ; end3 - start3