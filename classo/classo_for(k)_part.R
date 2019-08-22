k = 3
# threshold near-zero rhos to zero and stop algorithm
if(rho_path[k-1] <= 1e-4) {
  rho_path[k-1] = 0
  break
}

# calculate derivative for coefficients and multipliers
# construct matrix
activeCoeffs = which(setActive)
inactiveCoeffs = which(!setActive)
idxIneqBorder = which(setIneqBorder)

# 여기 계산 불안정 - 제약조건이 1개만(등호 or 부등호) 있을 때*****
M = cbind(H[activeCoeffs, activeCoeffs], Aeq[, activeCoeffs], 
          t(Aineq[setIneqBorder, activeCoeffs]))
M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs], Aineq[idxIneqBorder, activeCoeffs])

# calculate derivative
tryCatch(
  expr = {
    dir = dirsgn * solve(M, rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
  }, 
  finally = {
    dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
    # julia는 pinv 함수 씀
  }
)

if(any(is.nan(dir))) {
  dir = -(MASS::ginv(M) %*% rbind(subgrad[setActive], matrix(rep(0, neq + nIneqBorder), ncol = 1)))
}

# calculate the derivative for rho * subgradient
# 여기 계산 불안정 - 제약조건이 1개만 있을 때*****
temp = cbind(H[inactiveCoeffs, activeCoeffs], Aeq[, inactiveCoeffs])
if(nineq == 0) {
  dirSubgrad = - temp %*% dir
} else {
  dirSubgrad = - cbind(temp, Aineq[idxIneqBorder, inactiveCoeffs]) %*% dir  
}

### check additional events related to potential subgraient violations

## inactive coefficients moving too slowly
# negative subgradient
inactSlowNegIdx = which((1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                          subgrad[!setActive] <= (1 * dirsgn + 1e-8) &
                          1 * dirsgn < dirSubgrad)

# positive subgradient
inactSlowPosIdx = which((-1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                          subgrad[!setActive] <= (-1 * dirsgn + 1e-8) &
                          -1 * dirsgn > dirSubgrad)

## "active" coeficients estimated as 0 with potential sign mismatch #%
# positive subgrad but negative derivative
signMismatchPosIdx = which((0 - 1e-8) <= subgrad[setActive] &
                             subgrad[setActive] <= (1 + 1e-8) &
                             dirsgn * dir[1:nActive] <= (0 - 1e-8) &
                             beta_path[activeCoeffs, k-1] == 0)

# Negative subgradient but positive derivative
signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                             subgrad[setActive] <= (0 + 1e-8) &
                             dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                             beta_path[activeCoeffs, k-1] == 0)

# reset violation counter (to avoid infinite loops)
violateCounter = 0
