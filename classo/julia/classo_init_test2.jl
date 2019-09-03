# setting
cd("C:/Users/dpelt/OneDrive - 서울시립대학교/Desktop/Mayson/UOS_graduate/Constrained LASSO")
#cd("C:/Users/dpelt/AppData/Local/Julia-1.2.0\bin")
#cd("C:/Users/dpelt/AppData/Local/Julia-1.2.0/classo")
#using Pkg
#Pkg.add("Convex")
#Pkg.add("Statistics")
#Pkg.add("LinearAlgebra")
#Pkg.add("ECOS")
#Pkg.add("JuMP")
#Pkg.add("GLPk")
using Convex, Random, Statistics, LinearAlgebra, ECOS
#using classo_init

# basic setting ---------------------------------------------
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

#n = size(X, 1)
#p = size(X, 2)
#m = size(Aeq, 1)

# Initialization ---------------------------------------------
# path_init function part
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
rho_max = target_value[1]
lambda_max = target_value[(p+2):end]

# violation check
idx1 = findall(x -> x <= 1e-4, abs.((-X' * y + rho_max .* ones(p)) - Aeq' * lambda_max))
idx2 = findall(x -> x <= 1e-4, abs.((-X' * y - rho_max .* ones(p)) - Aeq' * lambda_max))
active_set = union(idx1, idx2)
num_active = size(active_set, 1)

# Iterate ---------------------------------------------
max_iter = 100
beta_path = zeros(p, max_iter)
rho_path = zeros(max_iter)
lambda_path = zeros(m, max_iter)

# save path
beta_path[:, 1] = zeros(p, 1)
rho_path[1] = rho_max
lambda_path[:, 1] = lambda_max

# find direction
H = X' * X
M = vcat(hcat(H[active_set, active_set], (Aeq[:, active_set])'),
    hcat(Aeq[:, active_set], zeros(m, m)))
resid = y - X * beta_path[:, 1]
subgrad = (- X' * resid - Aeq' * lambda_path[:, 1]) ./ rho_path[1]
sub_active_idx = findall(x -> 1 - abs(x) <= 1e-06, subgrad)
subgrad[sub_active_idx] = sign.(subgrad[sub_active_idx]) 
S = vcat(subgrad[sub_active_idx], zeros(m ,1))
b_m = pinv(M) * S
b = b_m[1:num_active ,1]
m = b_m[(num_active+1):end ,1]

# find delta_rho
delta_rho = 0
# c1
if delta_rho == 0
    println("finish!")
end
# c2

# c3 -> not using hitting the knot