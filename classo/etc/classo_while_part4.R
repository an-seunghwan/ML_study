# Monitor & fix condition 3 violations
while(length(signMismatchNegIdx) > 0) {
  ## Identify and move problem coefficient
  # indices corresponding to inactive coefficients
  activeCoeffs = which(setActive)
  # identify problem coefficient
  viol_coeff = activeCoeffs[signMismatchNegIdx]
  # put problem coefficient back into inactive set;
  setActive[viol_coeff] = F
  # determine new number of active coefficients
  nActive = sum(setActive != 0)
  # determine number of active/binding inequality constraints
  nIneqBorder = sum(setIneqBorder != 0)
  
  ## Recalculate derivative for coefficients & multipliers #%
  # construct matrix
  activeCoeffs = which(setActive)
  inactiveCoeffs = which(!setActive)
  idxIneqBorder = which(setIneqBorder)
  
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
  
  # calculate the derivative for rho * subgradient
  # 여기 계산 불안정 - 제약조건이 1개만 있을 때*****
  temp = cbind(H[inactiveCoeffs, activeCoeffs], Aeq[, inactiveCoeffs])
  if(nineq == 0) {
    dirSubgrad = - temp %*% dir
  } else {
    dirSubgrad = - cbind(temp, Aineq[idxIneqBorder, inactiveCoeffs]) %*% dir  
  }
  
  ## Misc. housekeeping #%
  # check for violations again
  # Negative subgradient but positive derivative
  signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                               subgrad[setActive] <= (0 + 1e-8) &
                               dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                               beta_path[activeCoeffs, k-1] == 0)
  
  # update violation counter
  violateCounter = violateCounter + 1
  # break loop if needed
  if(violateCounter >= maxiters) {
    break
  }
}