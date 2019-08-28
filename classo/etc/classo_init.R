rm(list=ls())
gc()
set.seed(520)
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")
# setting
n = 100
p = 10
m = 2

# data
X = matrix(rnorm(n*p), nrow = n)
X.m = apply(X, 2, mean)
X.m = matrix(X.m, nrow = n, ncol = p, byrow = T)
X = X - X.m
true_b = rep(1, p)
y = X %*% true_b + rnorm(n)
y = y - mean(y)

# constraints
Aeq = matrix(sample(seq(-1, 2, by = 1), m * p, replace = T), nrow = m)
# Aeq = matrix(c(1,1,0,
               # -1,0,1,
               # -1,-1,1), nrow = m, byrow = T)
beq = matrix(rep(0, m), nrow = m)

# setting
zero_p = matrix(rep(0, p), ncol = 1)
zero_m = matrix(rep(0, m), ncol = 1)
zero_pm = matrix(rep(0, p*m), nrow = p)
zero_pp = matrix(rep(0, p*p), nrow = p)
one_p = matrix(rep(1, p), ncol = 1)
one_m = matrix(rep(1, m), ncol = 1)
I_p = matrix(diag(rep(1, p)), nrow = p)

# design matrix of Equality
Aeq1 = matrix(c(1, t(zero_p), t(zero_m)), nrow = 1)
Aeq2 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(zero_p, I_p, zero_pm),
             c(0, t(zero_p), t(zero_m)))
Aeq3 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(zero_p, zero_pp, t(Aeq)),
             c(0, t(zero_p), t(zero_m)))
Aeq4 = rbind(0, -t(X)%*%y, 0)
Aeq5 = rbind(c(0, t(zero_p), t(zero_m)),
             cbind(one_p, zero_pp, zero_pm),
             c(0, t(zero_p), t(zero_m)))

# solve
target = Variable(1 + p + m)
constraints = list(Aeq2 %*% target == Aeq3 %*% target,
                   Aeq2 %*% target <= Aeq4 + Aeq5 %*% target,
                   Aeq2 %*% target >= Aeq4 - Aeq5 %*% target,
                   Aeq1 %*% target >= 0)
objective = Minimize(Aeq1 %*% target)
problem = Problem(objective, constraints)
result = solve(problem)

# result
target = result$getValue(target)
rho_min = target[1, ]
z = target[2:(1+p), ,drop=F]
lambda = target[(p+2):nrow(target), ,drop=F]
rho_min; z; lambda
Aeq

# check
sum(abs(z - t(Aeq) %*% lambda)) < 1e-6

# violation check : rho is decreasing
t(Aeq) %*% lambda
((-t(X) %*% y + rho_min) - (t(Aeq) %*% lambda)) < 1e-8
which.min((-t(X) %*% y + rho_min) - (t(Aeq) %*% lambda))

# maximum of lowerbound
-t(X) %*% y + rho_min
-t(X) %*% y - rho_min
min(-t(X) %*% y + rho_min)
max(-t(X) %*% y - rho_min)

t(Aeq) %*% lambda

t(Aeq[, 1]) %*% lambda

# init activeset?
max.j = which.max(-t(X) %*% y - rho_min) # maximum of lowerbound
min.j = which.min(-t(X) %*% y + rho_min) # minimum of upperbound

# plot
## 이거 3차원으로 (l1, l2, affine value) 그려서 -t(X) %*% y - rho_min의 scalar 범위 확인해서 그려봐야 할듯!
l1 = seq(-300, 100, length.out = 10000)
l2 = seq(-300, 100, length.out = 10000)

plot(l1, l2, ylab = "", xlab = "", type = "n")
abline(v = (-t(X) %*% y - rho_min))
abline(v = (-t(X) %*% y + rho_min))
abline(h = (-t(X) %*% y - rho_min))
abline(h = (-t(X) %*% y + rho_min))
for(i in 1:p) lines(l1, t(Aeq[, i])[1] * l1 + t(Aeq[, i])[2] * l2, ylab = "", xlab = "", col = i, lwd = 1)

plot(l1, l2, ylab = "", xlab = "", type = "n")
abline(v = max(-t(X) %*% y - rho_min), h = max(-t(X) %*% y - rho_min))
abline(v = min(-t(X) %*% y + rho_min), h = min(-t(X) %*% y + rho_min))
for(i in 1:p) lines(l1, t(Aeq[, i])[1] * l1 + t(Aeq[, i])[2] * l2, ylab = "", xlab = "", col = i, lwd = 1)

plot(l1, l2, ylab = "", xlab = "", type = "n")
abline(v = max(-t(X) %*% y - rho_min), h = max(-t(X) %*% y - rho_min))
abline(v = min(-t(X) %*% y + rho_min), h = min(-t(X) %*% y + rho_min))
for(i in c(max.j, min.j)) lines(l1, t(Aeq[, i])[1] * l1 + t(Aeq[, i])[2] * l2, ylab = "", xlab = "", col = i, lwd = 1)

plot(l1, l2, ylab = "", xlab = "", type = "n")
for(i in 1:p) {
  abline(v = (-t(X) %*% y - rho_min)[i, ], h = (-t(X) %*% y - rho_min)[i, ], col = i)
  abline(v = (-t(X) %*% y + rho_min)[i, ], h = (-t(X) %*% y + rho_min)[i, ], col = i)
  lines(l1, t(Aeq[, i])[1] * l1 + t(Aeq[, i])[2] * l2, ylab = "", xlab = "", col = i, lwd = 1)
}

plot(l1, t(Aeq[, 1])[1] * l1 + t(Aeq[, 1])[2] * l2, ylab = "", xlab = "", col = 1, lwd = 1, type = 'l')
for(i in 2:p) lines(l1, t(Aeq[, i])[1] * l1 + t(Aeq[, i])[2] * l2, ylab = "", xlab = "", col = i, lwd = 1)
for(i in 2:p) lines(l1, t(Aeq[, i])[1] * l1 + t(Aeq[, i])[2] * l2, ylab = "", xlab = "", col = i, lwd = 1)









