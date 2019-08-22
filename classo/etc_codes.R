
##########################################################
c = matrix(rep(1, p), nrow = 1)
ECOS_csolve(c = c, G = Aineq, h = matrix(bineq, nrow = p), dims = list(l = dim(X)[1] * 2, q = NULL, e = 0), A = Aeq, b = beq)


#############################################################
x = Variable(p)

p = JuliaObject(p)
Aeq = JuliaObject(Aeq)
beq = JuliaObject(beq)
Aineq = JuliaObject(Aineq)
bineq = JuliaObject(bineq)

problem = minimize(sum(abs(x)), Aeq %*% x == beq, Aineq %*% x <= bineq)

julia_assign("x", Variable(p))
julia_assign("p", p)
julia_assign("Aeq", Aeq)
julia_assign("beq", beq)
julia_assign("Aineq", Aineq)
julia_assign("bineq", bineq)

julia_command("problem = minimize(dot(ones(p), abs(x)), Aeq * x == beq, Aineq * x <= bineq)")

julia_call("minimize", problem)



problem = julia_call("minimize", sum(abs(x))
                     , Aeq * x == beq, Aineq * x <= bineq)



# (앞과 동일한 식 아닌가? 굳이 2번?)
# problem = minimize(dot(ones(eltype(X), size(X, 2)), abs(x)))
# if !isempty(Aeq)
# problem.constraints += Aeq * x == beq
# end
# if !isempty(Aineq)
# problem.constraints += Aineq * x <= bineq
# end

# TT = STDOUT # save original STDOUT stream
# redirect_stdout()
a =
# redirect_stdout(TT) # restore STDOUT


  #############################################################################################
p = dim(X)[2]

x = Variable(p)
objective = Minimize(sum(abs(x)))
constraints = list(Aineq %*% x <= bineq, Aeq %*% x == beq) # 제약조건을 반드시 부등호 -> 등호 조건 순서대로 기입
problem = Problem(objective, constraints)
result = solve(problem)

beta = result$getValue(x)

lambda_eq = matrix(0, 0, 1)
mu_ineq = matrix(0, 0, 1)

for(i in 1:min(2, length(problem@constraints))) {
  if(canonicalize(problem@constraints[[i]])[[2]][[1]]$class == "LinEqConstr") {
    lambda_eq = result$getDualValue(problem@constraints[[i]])
  } else if(canonicalize(problem@constraints[[i]])[[2]][[1]]$class == "LinLeqConstr") {
    mu_ineq = result$getDualValue(problem@constraints[[i]])
  }
}
# Error in subset.default(x, subset) : 'subset' must be logical (무슨 오류?)

# for(i in 1:min(2, length(constraints))) {
#   if(str_detect(as.character(constraints[[i]]), "==")) {
#     lambda_eq = result$getDualValue(constraints[[i]])
#   } else if(str_detect(as.character(constraints[[i]]), "<=")) {
#     mu_ineq = result$getDualValue(constraints[[i]])
#   }
# }

if(sum(!is.empty(mu_ineq))) {
  mu_ineq[mu_ineq < 0] = 0
}

setActive = abs(beta) > 1e-4 | !penidx
beta[!setActive] = 0

resid = y - X %*% beta
subgrad = t(X) %*% resid - t(Aeq) %*% lambda_eq - t(Aineq) %*% mu_ineq
lo_max = max(abs(subgrad))
ind_lo_max = which.max(abs(subgrad))

###################################################################################################

#################################################################
n = dim(X)[1]
p = dim(X)[2]

# function 인자값 setting
warmstart = F
rho = rho_path[1]
# penwt = penidx
beta_hat = matrix(rep(0, p*length(rho)), nrow = p, ncol = length(rho))
optval_vec = rep(0, length(rho))
prob_vec = list()

beta = Variable(p)
loss = (1 / 2) * sum_squares(sqrt(obswt) * (y - X %*% beta)) # rhoss term
pen = penwt %*% abs(beta)

# get constraints
if(neq > 0 & nineq > 0) {
  constraints = list(Aeq %*% x == beq, Aineq %*% x <= bineq)
} else if(neq > 0 & nineq == 0) {
  constraints = list(Aeq %*% x == beq)
} else if(neq == 0 & nieq > 0) {
  constraints = list(Aineq %*% x <= bineq)
}

for(i in 1:length(rho)) {
  rho_i = rho[i]
  if(rho_i == Inf) {
    objective = Minimize(pen)
    problem = Problem(objective, constraints)
  } else if(rho_i <= 0) {
    objective = Minimize(loss)
    problem = Problem(objective, constraints)
  } else {
    objective = Minimize(loss + rho_i * pen)
    problem = Problem(objective, constraints)
  }
  
  if(warmstart) {
    result = solve(problem, warm_start = i) # warm_start 잘 모르겠음...
    prob_vec[[i]] = problem
    optval_vec[i] = result$value
  } else {
    result = solve(problem) 
    prob_vec[[i]] = problem
    optval_vec[i] = result$value
  }
  
  if(length(rho) == 1) {
    return(list(beta_value = result$getValue(beta),
                problem.optval = result$value,
                problem = problem))
  }
  
  beta_hat[, i] = result$getValue(x)
}

return(list(beta_hat = beta_hat,
            optval_vec = optval_vec,
            prob_vec = prob_vec))
