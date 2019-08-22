rm(list = ls())
gc()
# if(!require("convexjlr")) install.packages("convexjlr")
# library("convexjlr")
# if(!require("ECOSolveR")) install.packages("ECOSolveR")
# library("ECOSolveR")
# if(!require("cccp")) install.packages("cccp")
# library("cccp")
# if(!require("lpSolve")) install.packages("lpSolve")
# library("lpSolve")
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")
if(!require("rapportools")) install.packages("rapportools")
library("rapportools")
# if(!require("stringr")) install.packages("stringr")
# library("stringr")
# julia <- julia_setup()
# convex_setup()

### set up
n = 200
p = 20

### truth with sum constrant sum(beta) = 0
# β = zeros(p)
# β[1:round(Int, p / 4)] = 0
# β[(round(Int, p / 4) + 1):round(Int, p / 2)] = 1
# β[(round(Int, p / 2) + 1):round(Int, 3p / 4)] = 0
# β[(round(Int, 3p / 4) + 1):p] = -1

### generate data
X=as.matrix(read.csv("C:/Julia/classo/src/X.csv", header = F)) # n * p
y=as.matrix(read.csv("C:/Julia/classo/src/y.csv", header = F)) # n

### equality constraints
Aeq = matrix(rep(1, p), nrow = 1) # 1×20 Array{Float64,2}
beq = matrix(0, nrow = 1) # 1-element Array{Float64,1}
penwt = rep(1, p) # 20-element Array{Float64,1}

# inequality constraints
Aineq = -diag(rep(1, p))
bineq = rep(0, p)

# globla variable로 처리하면 좋을 것 같은 변수 *****
# : Aeq, beq, Aineq, bineq, ρridge, penidx

### lsq_classopath 내부
# X      :: AbstractMatrix{T},
# y      :: AbstractArray{T};
# Aeq    :: AbstractMatrix = zeros(T, 0, size(X, 2)),
# beq    :: Union{AbstractArray, Number} = zeros(T, size(Aeq, 1)),
# Aineq  :: AbstractMatrix = zeros(T, 0, size(X, 2)),
# bineq  :: Union{AbstractArray, Number} = zeros(T, size(Aineq, 1)),
# ρridge :: Number = zero(T),
# penidx :: Array{Bool} = fill(true, size(X, 2)),
# solver = ECOSSolver(maxit=10e8, verbose=0)

n = dim(X)[1]
p = dim(X)[2]

lo_ridge = 0
penidx = rep(T, p)


if(n < p) {
  print("Adding a small ridge penalty (default is 1e-4) since n < p")
  if(lo_ridge <= 0) {
    print("ρridge must be positive, switching to default value (1e-4)")
    lo_ridge = 1e-4
  }
  # create augmented data
  y = rbind(y, matrix(rep(0, p), nrow = p))
  X = rbind(X, sqrt(lo_ridge) * diag(rep(1, p)))
  # record original number of observations
  # n_orig = n
} else {
  # make sure X is full column rank
  qrfact_X = qr(X)
  R = qr.R(qrfact_X)
  rankX = sum(abs(diag(R)) > (abs(R[1,1]) * max(n,p) * (.Machine$double.eps) ^ 2))  # 4.930380657631324e-32 *****
  
  if(rankX != p) {
    print("Adding a small ridge penalty (default is 1e-4) since X is rank deficient")
    if(lo_ridge <= 0) {
      print("ρridge must be positive, switching to default value (1e-4)")
      lo_ridge = 1e-4
    }
    # create augmented data
    y = rbind(y, matrix(rep(0, p), nrow = p))
    X = rbind(X, sqrt(lo_ridge) * diag(rep(1, p)))
  }
}

# allocate variables along path
neq = dim(Aeq)[1]
nineq = dim(Aineq)[1]
maxiters = 5 * (p + nineq) # max number of path segments to consider
beta_path = matrix(rep(0, p * maxiters), nrow = p) 
lambda_path = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
mu_path = matrix(rep(0, nineq * maxiters), nrow = nineq) # dual variables for inequality
lo_path = rep(0, maxiters) # tuning parameter
df_path = rep(Inf, maxiters) # degree of freedom
objval_path = rep(0, maxiters) # objective value
violation_path = rep(Inf, maxiters) 

### initialization
# use LP to find ρmax
H = t(X) %*% X 
#sense = [repmat(['='], neq); repmat(['<'], nineq)]
#β, _, problem = lsq_constrsparsereg(X, y, Inf;      # why necessary?
#          [Aeq; Aineq], sense, [beq; bineq], lb, ub)#

l = find_lo_max()
lo_path[1] = l$lo_max
idx = l$ind_lo_max
























