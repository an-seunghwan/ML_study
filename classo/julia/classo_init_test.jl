# setting
#cd("C:/Users/dpelt/OneDrive - 서울시립대학교/Desktop/Mayson/UOS_graduate/Constrained LASSO")
#cd("C:\Users\dpelt\AppData\Local\Julia-1.2.0\bin")
using Pkg
#Pkg.add("Convex")
#Pkg.add("Statistics")
#Pkg.add("LinearAlgebra")
#Pkg.add("ECOS")
#Pkg.add("JuMP")
Pkg.add("GLPk")
using Convex, Random, Statistics, LinearAlgebra, JuMP, GLPK

# parameters
Random.seed!(520)
n, p, m = 100, 10, 5
X = randn(n, p)
X_mean = mean(X, dims = 1)
X = X - repeat(X_mean, n)
true_b = [1.0,2.0,3.0,2.0,1.0,3.0,2.0,1.0,2.0,1.0]
y = X * true_b + randn(n)

# constraints
Aeq = randn(m, p)
beq = zeros(m)

function path_init(X, y, Aeq, beq) 
    n = size(X, 1)
    p = size(X, 2)
    m = size(Aeq, 1)

    D1 = hcat(1, zeros(p)', zeros(m)')
    D2 = vcat(hcat(0, zeros(p)', zeros(m)'), 
        hcat(zeros(p), Matrix(I, p, p), zeros(p, m)), 
        hcat(0, zeros(p)', zeros(m)'))
    D3 = vcat(hcat(0, zeros(p)', zeros(m)'), 
        hcat(zeros(p), zeros(p, p), Aeq'), 
        hcat(0, zeros(p)', zeros(m)'))
    D4 = vcat(0, -X' * y, 0)
    D5 = vcat(hcat(0, zeros(p)', zeros(m)'), 
        hcat(ones(p), zeros(p, p), zeros(p, m)), 
        hcat(0, zeros(p)', zeros(m)'))

    # solve
    target = Variable(1 + p + m)
    problem = minimize(D1 * target,
        D2 * target == D3 * target,
        D2 * target <= D4 + D5 * target,
        D2 * target >= D4 - D5 * target,
        D1 * target >= 0)
    solver = ECOSSolver(maxit=10e8, verbose=0)
    solve!(problem, solver)
    target_value = target.value
    rho_max = target_value[1, :]
    lambda_max = target_value[(p+2):end, :]

    # violation check
    idx1 = findall(abs.((-X' * y + rho_max * ones(p)) - Aeq' * lambda_max) .<= 1e-4)
    idx2 = findall(abs.((-X' * y - rho_max * ones(p)) - Aeq' * lambda_max) .<= 1e-4)
    active_set = union(idx1, idx2)

    # Aeq rank 확인하는 부분 추가하면 좋을듯!

    return rho_max, lambda_max, active_set
end

a, b, c = path_init(X, y, Aeq, beq)

n = size(X, 1)
p = size(X, 2)
m = size(Aeq, 1)

D1 = hcat(1, zeros(p)', zeros(m)')
D2 = vcat(hcat(0, zeros(p)', zeros(m)'), 
    hcat(zeros(p), Matrix(I, p, p), zeros(p, m)), 
    hcat(0, zeros(p)', zeros(m)'))
D3 = vcat(hcat(0, zeros(p)', zeros(m)'), 
    hcat(zeros(p), zeros(p, p), Aeq'), 
    hcat(0, zeros(p)', zeros(m)'))
D4 = vcat(0, -X' * y, 0)
D5 = vcat(hcat(0, zeros(p)', zeros(m)'), 
    hcat(ones(p), zeros(p, p), zeros(p, m)), 
    hcat(0, zeros(p)', zeros(m)'))

# solve
target = Variable(1 + p + m)
problem = minimize(D1 * target,
    D2 * target == D3 * target,
    D2 * target <= D4 + D5 * target,
    D2 * target >= D4 - D5 * target,
    D1 * target >= 0)
solver = ECOSSolver(maxit=10e8, verbose=0)
solve!(problem, solver)
target_value = target.value
rho_max = target_value[1, :]
lambda_max = target_value[(p+2):end, :]

# violation check
idx1 = findall(abs.((-X' * y + rho_max * ones(p)) - Aeq' * lambda_max) .<= 1e-4)
idx2 = findall(abs.((-X' * y - rho_max * ones(p)) - Aeq' * lambda_max) .<= 1e-4)
active_set = union(idx1, idx2)