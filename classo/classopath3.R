classopath_modified = function(X, y, Aeq, beq, Aineq, bineq,
                      penwt = rep(1, p),
                      rho_ridge = 0,
                      penidx = rep(T, p),
                      choose_one = F) {
  n = dim(X)[1]
  p = dim(X)[2]
  
  if(n < p) {
    print("Adding a small ridge penalty (default is 1e-4) since n < p")
    if(rho_ridge <= 0) {
      print("ρridge must be positive, switching to default value (1e-4)")
      rho_ridge = 1e-4
    }
    # create augmented data
    y = rbind(y, matrix(rep(0, p), nrow = p))
    X = rbind(X, sqrt(rho_ridge) * diag(rep(1, p)))
    # record original number of observations
    # n_orig = n
  } else {
    # make sure X is full column rank
    qrfact_X = qr(X)
    R = qr.R(qrfact_X)
    rankX = sum(abs(diag(R)) > (abs(R[1,1]) * max(n,p) * (.Machine$double.eps) ^ 2))  # 4.930380657631324e-32 *****
    
    if(rankX != p) {
      print("Adding a small ridge penalty (default is 1e-4) since X is rank deficient")
      if(rho_ridge <= 0) {
        print("ρridge must be positive, switching to default value (1e-4)")
        rho_ridge = 1e-4
      }
      # create augmented data
      y = rbind(y, matrix(rep(0, p), nrow = p))
      X = rbind(X, sqrt(rho_ridge) * diag(rep(1, p)))
    }
  }
  
  # alrhocate variables arhong path
  neq = dim(Aeq)[1]
  nineq = dim(Aineq)[1]
  maxiters = 5 * (p + nineq) # max number of path segments to consider
  beta_path = matrix(rep(0, p * maxiters), nrow = p) 
  lambda_patheq = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
  mu_pathineq = matrix(rep(0, nineq * maxiters), nrow = nineq) # dual variables for inequality
  rho_path = rep(0, maxiters) # tuning parameter
  df_path = rep(Inf, maxiters) # degree of freedom
  objval_path = rep(0, maxiters) # objective value
  violation_path = rep(Inf, maxiters) 
  
  ### initialization
  H = t(X) %*% X 
  # find the maximum ρ (starting value)
  l = init_path(X, y, Aeq, beq) # no inequality constraints
  rho_path[1] = l$rho_max
  lambda_patheq[, 1] = l$lambda_max
  if(!l$flag) return(0) # check lambda is unique
  setActive = matrix(rep(F, p), ncol = 1)
  setActive[l$activeset, ] = T
  
  #
  if(choose_one == T) {
    # choose one predictor most correlated
    setActive = matrix(rep(F, p), ncol = 1)
    idx = l$activeset[which.max(t(X[, l$activeset]) %*% y)]
    setActive[idx, ] = T
  }
  
  beta_path[, 1] = 0
  objval_path[1] = sum(t(X) %*% y) / 2 # beta is zero
  
  #
  mu_pathineq[mu_pathineq < 0] = 0
  residIneq = Aineq %*% beta_path[, 1] - bineq
  setIneqBorder = residIneq == 0
  nIneqBorder = sum(setIneqBorder != 0)
  
  # initialize subgradient vector (stationarity condition)
  resid = y - X %*% beta_path[, 1]
  if(neq > 0 & nineq > 0) {
    subgrad = - t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1] - t(Aineq) %*% mu_pathineq[, 1]
  } else if(neq > 0 & nineq == 0) {
    subgrad = - t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]
  } else if(neq == 0 & nineq > 0) {
    subgrad = - t(X) %*% resid - t(Aineq) %*% mu_pathineq[, 1]
  }
  
  # subgradient of beta
  # subgrad[setActive] = sign(beta_path[setActive, 1])
  # subgrad[!setActive] = subgrad[!setActive] / rho_path[1]
  subgrad = subgrad / rho_path[1] # 처음에는 beta가 다 0이니까 sign말고 직접 계산?
  nActive = sum(setActive != 0)
  
  # calculate degrees of freedom
  # rankAeq = Matrix::rankMatrix(Aeq) ### Brian's comment: need to make it more efficient
  rankAeq = ifelse(dim(Aeq)[1] > 0, Matrix::rankMatrix(Aeq), 0)
  df_path[1] = nActive - rankAeq - nIneqBorder
  
  # set initial violations counter to 0
  violation_path[1] = 0
  
  # sign for path direction (originally went both ways, but increasing was retired -> only decreasing)
  dirsgn = -1
  
  ####################################
  ### main loop for path following ###
  ####################################
  
  for(k in 2:maxiters) {
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
    
    M = cbind(H[activeCoeffs, activeCoeffs, drop=F], t(Aeq[, activeCoeffs, drop=F]), 
              t(Aineq[setIneqBorder, activeCoeffs, drop=F]))
    M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
    M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs, drop=F], 
                                                                    Aineq[idxIneqBorder, activeCoeffs, drop=F])
    
    # calculate derivative (of beta, lambda, mu)
    if(min(eigen(M)$values) > 0) {
      dir = -dirsgn * solve(M, rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
    } else {
      dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
    }
    
    # calculate the derivative for rho * subgradient
    dirSubgrad = - cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]),
                         t(Aineq[idxIneqBorder, inactiveCoeffs, drop=F])) %*% dir
    
    ### check additional events related to potential subgraient violations
    
    ## inactive coefficients moving too slowly
    # negative subgradient
    inactSlowNegIdx = which((1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                              subgrad[!setActive] <= (1 * dirsgn + 1e-8) & # 매우 좁은 범위의 값으로 -1과 값 확인
                              1 * dirsgn < dirSubgrad)
    
    # positive subgradient
    inactSlowPosIdx = which((-1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                              subgrad[!setActive] <= (-1 * dirsgn + 1e-8) &
                              -1 * dirsgn > dirSubgrad)
    
    ## "active" coeficients estimated as 0 with "potential" sign mismatch #% ???
    # positive subgrad but negative derivative
    signMismatchPosIdx = which((0 - 1e-8) <= subgrad[setActive] &
                                 subgrad[setActive] <= (1 + 1e-8) & # 처음 2 조건: beta가 양수임
                                 dirsgn * dir[1:nActive] <= (0 - 1e-8) & # 변화량: 음수
                                 beta_path[activeCoeffs, k-1] == 0) # 이전 beta: 0 
    
    # Negative subgradient but positive derivative
    signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                                 subgrad[setActive] <= (0 + 1e-8) &
                                 dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                                 beta_path[activeCoeffs, k-1] == 0)
    
    # reset violation counter (to avoid infinite loops)
    violateCounter = 0
    
    # outer while loop for checking all conditions together
    while(length(inactSlowNegIdx) > 0 | length(inactSlowPosIdx) > 0 |
          length(signMismatchPosIdx) > 0 | length(signMismatchNegIdx) > 0) {
      
      # monitor and fix condition 1 violations
      while(length(inactSlowNegIdx) > 0) {
        ## Identify and move problem coefficient
        # indices corresponding to inactive coefficients
        inactiveCoeffs = which(!setActive)
        # identify problem coefficient
        viol_coeff = inactiveCoeffs[inactSlowNegIdx]
        setActive[viol_coeff] = T
        
        nActive = sum(setActive != 0)
        nIneqBorder = sum(setIneqBorder != 0)
        
        activeCoeffs = which(setActive)
        inactiveCoeffs = which(!setActive)
        idxIneqBorder = which(setIneqBorder)
        
        M = cbind(H[activeCoeffs, activeCoeffs, drop=F], t(Aeq[, activeCoeffs, drop=F]), 
                  t(Aineq[setIneqBorder, activeCoeffs, drop=F]))
        M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
        M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs, drop=F], 
                                                                        Aineq[idxIneqBorder, activeCoeffs, drop=F])
        
        if(min(eigen(M)$values) > 0) {
          dir = -dirsgn * solve(M, rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        } else {
          dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        }
        
        # calculate the derivative for rho * subgradient
        dirSubgrad = - cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]),
                             t(Aineq[idxIneqBorder, inactiveCoeffs, drop=F])) %*% dir
        
        ## Misc. housekeeping #%
        # check for violations again
        
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
        
        # update violation counter
        violateCounter = violateCounter + 1
        # break loop if needed
        if(violateCounter >= maxiters) {
          break
        }
      }
      
      # Monitor & fix subgradient condition 2 violations
      while(length(inactSlowPosIdx) > 0) {
        ## Identify and move problem coefficient
        # indices corresponding to inactive coefficients
        inactiveCoeffs = which(!setActive)
        # identify problem coefficient
        viol_coeff = inactiveCoeffs[inactSlowPosIdx]
        setActive[viol_coeff] = T
        
        nActive = sum(setActive != 0)
        nIneqBorder = sum(setIneqBorder != 0)
        
        activeCoeffs = which(setActive)
        inactiveCoeffs = which(!setActive)
        idxIneqBorder = which(setIneqBorder)
        
        M = cbind(H[activeCoeffs, activeCoeffs, drop=F], t(Aeq[, activeCoeffs, drop=F]), 
                  t(Aineq[setIneqBorder, activeCoeffs, drop=F]))
        M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
        M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs, drop=F], 
                                                                        Aineq[idxIneqBorder, activeCoeffs, drop=F])
        
        # calculate derivative
        if(min(eigen(M)$values) > 0) {
          dir = -dirsgn * solve(M, rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        } else {
          dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        }
        
        # calculate the derivative for rho * subgradient
        dirSubgrad = - cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]),
                             t(Aineq[idxIneqBorder, inactiveCoeffs, drop=F])) %*% dir
        
        ## Misc. housekeeping #%
        # check for violations again
        
        ## inactive coefficients moving too slowly
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
        
        # update violation counter
        violateCounter = violateCounter + 1
        # break loop if needed
        if(violateCounter >= maxiters) {
          break
        }
      }
      
      # Monitor & fix condition 3 violations
      while(length(signMismatchPosIdx) > 0) {
        ## Identify and move problem coefficient
        # indices corresponding to inactive coefficients
        activeCoeffs = which(setActive)
        # identify problem coefficient
        viol_coeff = activeCoeffs[signMismatchPosIdx]
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
        
        M = cbind(H[activeCoeffs, activeCoeffs, drop=F], t(Aeq[, activeCoeffs, drop=F]), 
                  t(Aineq[setIneqBorder, activeCoeffs, drop=F]))
        M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
        M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs, drop=F], 
                                                                        Aineq[idxIneqBorder, activeCoeffs, drop=F])
        
        # calculate derivative
        if(min(eigen(M)$values) > 0) {
          dir = -dirsgn * solve(M, rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        } else {
          dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        }
        
        # calculate the derivative for rho * subgradient
        dirSubgrad = - cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]),
                             t(Aineq[idxIneqBorder, inactiveCoeffs, drop=F])) %*% dir
        
        ## Misc. housekeeping #%
        # check for violations again
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
        
        # update violation counter
        violateCounter = violateCounter + 1
        # break loop if needed
        if(violateCounter >= maxiters) {
          break
        }
      }
      
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
        
        M = cbind(H[activeCoeffs, activeCoeffs, drop=F], t(Aeq[, activeCoeffs, drop=F]), 
                  t(Aineq[setIneqBorder, activeCoeffs, drop=F]))
        M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
        M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs, drop=F], 
                                                                        Aineq[idxIneqBorder, activeCoeffs, drop=F])
        
        # calculate derivative
        if(min(eigen(M)$values) > 0) {
          dir = -dirsgn * solve(M, rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        } else {
          dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        }
        
        # calculate the derivative for rho * subgradient
        dirSubgrad = - cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]),
                             t(Aineq[idxIneqBorder, inactiveCoeffs, drop=F])) %*% dir
        
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
      
      ## update violation trackers to see if any issues persist ##%
      activeCoeffs = which(setActive)
      inactiveCoeffs = which(!setActive)
      idxIneqBorder = which(setIneqBorder)
      
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
      
      # update violation counter
      violateCounter = violateCounter + 1
      # break loop if needed
      if(violateCounter >= maxiters) {
        break
      }
    }
    
    ###### after all while loop
    
    # store number of violations
    violation_path[k] = violateCounter
    
    # calculate derivative for residual inequality
    activeCoeffs = which(setActive)
    inactiveCoeffs = which(!setActive)
    idxIneqBorder = which(setIneqBorder)
    
    dirResidIneq = Aineq[which(!setIneqBorder), activeCoeffs, drop=F] %*% dir[1:nActive] 
    
    ### Determine rho for next event (via delta rho) ###% ?????
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
    if(length(idx) > 0) {
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
  }
  
  beta_path = beta_path[, 1:k-1]
  rho_path = rho_path[1:k-1]
  objval_path = objval_path[1:k-1]
  if(neq > 0) lambda_patheq = lambda_patheq[, 1:k-1]
  if(nineq > 0) mu_pathineq = mu_pathineq[, 1:k-1]
  df_path = df_path[1:k-1]
  df_path[df_path < 0] = 0
  violation_path = violation_path[k-1]
  
  return(list(beta_path = beta_path,
              rho_path = rho_path,
              objval_path = objval_path,
              lambda_patheq = lambda_patheq,
              mu_pathineq = mu_pathineq,
              df_path = df_path,
              violation_path = violation_path,
              steps = k-1))
}