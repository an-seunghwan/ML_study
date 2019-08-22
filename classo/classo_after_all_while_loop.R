###### after all while loop

# store number of violations
violation_path[k] = violateCounter

# calculate derivative for residual inequality
activeCoeffs = which(setActive)
inactiveCoeffs = which(!setActive)
idxIneqBorder = which(setIneqBorder)

dirResidIneq = Aineq[which(!setIneqBorder), activeCoeffs] %*% dir[1:nActive] # 계산 불안정

### Determine rho for next event (via delta rho) ###%
## Events based on coefficients changing activation status ##%
next_rho_beta = matrix(rep(Inf, p), ncol = 1)

# Active coefficient going inactive #%
next_rho_beta[setActive] = -dirsgn * beta_path[activeCoeffs, k-1] / dir[1:nActive]

# Inactive coefficient becoming positive #%
t1 = dirsgn * rho_path[k-1] * (1 - subgrad[inactiveCoeffs]) / (dirSubgrad - 1)
# threshold values hitting ceiling
t1[t1 < 1e-8] = Inf

# Inactive coefficient becoming negative #%
t2 = -dirsgn * rho_path[k-1] * (1 + subgrad[!setActive]) / (dirSubgrad + 1)
# threshold values hitting ceiling
t2[t2 < 1e-8] = Inf

# choose smaller delta rho out of t1 and t2
next_rho_beta[!setActive] = pmin(t1, t2) # elementwise min
next_rho_beta[(next_rho_beta <= 1e-8) | !penidx] = Inf

## Events based inequality constraints ##%
# clear previous values
next_rho_Ineq = matrix(rep(Inf, nineq), ncol = 1)

# Inactive inequality constraint becoming active #%
next_rho_Ineq[!setIneqBorder] = matrix(-dirsgn * residIneq[!setIneqBorder], sum(!setIneqBorder != 0), 1) /
  matrix(dirResidIneq, sum(!setIneqBorder != 0), 1)

# Active inequality constraint becoming deactive #%
if(length(mu_pathineq) > 0) {
  next_rho_Ineq[setIneqBorder] = -dirsgn * mu_pathineq[idxIneqBorder, k-1] / 
    matrix(dir[nActive + neq + 1, ], nIneqBorder, 1)
}

next_rho_Ineq[next_rho_Ineq <= 1e-8] = Inf

## determine next rho ##
# find smallest rho
chg_rho = min(rbind(next_rho_beta, next_rho_Ineq), na.rm = T)
# find all indices corresponding to this chgρ
idx = which(rbind(next_rho_beta, next_rho_Ineq) - chg_rho <= 1e-8)

# terminate path following if no new event found
if(chg_rho == Inf) {
  chg_rho = rho_path[k-1]
}

## Update values at new rho ##%
# move to next rho #%
# make sure next rho isn't negative
if(rho_path[k-1] + dirsgn * chg_rho < 0) {
  chg_rho = rho_path[k-1]
}

# calculate new value of rho
rho_path[k] = rho_path[k-1] + dirsgn * chg_rho

## Update parameter and subgradient values #%
# new coefficient estimates
activeCoeffs = which(setActive)

beta_path[activeCoeffs, k] = beta_path[activeCoeffs, k-1] + dirsgn * chg_rho * dir[1:nActive]
# force near-zero coefficients to be zero (helps with numerical issues)
beta_path[abs(beta_path[, k]) < 1e-12, k] = 0

# new subgradient estimates
subgrad[!setActive] = (rho_path[k-1] * subgrad[!setActive] + 
                         dirsgn * chg_rho * dirSubgrad) / rho_path[k]

# Update dual variables #%
# update lambda (lagrange multipliers for equality constraints)
if(length(lambda_patheq) > 0) {
  lambda_patheq[, k] = lambda_patheq[, k-1] + 
    dirsgn * chg_rho * matrix(dir[(nActive + 1):(nActive + neq), ], neq, 1)
}
# update mu (lagrange multipliers for inequality constraints)
if(length(mu_pathineq) > 0) {
  mu_pathineq[idxIneqBorder, k] = mu_pathineq[idxIneqBorder, k-1] + 
    dirsgn * chg_rho * matrix(dir[(nActive + neq + 1) : ncol(dir)], nIneqBorder, 1)
}
# update residual inequality
residIneq = Aineq * beta_path[, k] - bineq

## update sets ##%
for(j in 1:length(idx)) {
  curidx = idx[j]
  if(curidx <= p & setActive[curidx]) {
    # an active coefficient hits 0, or
    setActive[curidx] = F
  } else if(curidx <= p & !setActive[curidx]) {
    # a zero coefficient becomes nonzero
    setActive[curidx] = T
  } else if(curidx > p) {
    # an ineq on boundary becomes strict, or
    # a strict ineq hits boundary
    setIneqBorder[curidx - p] = !setIneqBorder[curidx - p]
  }
}

# determine new number of active coefficients
nActive = sum(setActive != 0)
# determine number of active/binding inequality constraints
nIneqBorder = sum(nIneqBorder != 0)

## Calcuate and store values of interest along the path #%
# calculate value of objective function
objval_path[k] = norm(y - X %*% beta_path[, k], type = "F")^2 / 2 + 
  rho_path[k] * sum(abs(beta_path[, k]))

# calculate degrees of freedom
df_path[k] = nActive - rankAeq - nIneqBorder
# break algorithm when df are exhausted
if(df_path[k] >= n) {
  break
}

# end of big for loop



































