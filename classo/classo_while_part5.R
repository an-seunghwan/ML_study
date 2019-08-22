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