# SVM #

using LinearAlgebra
using Random
using Printf

# read data
path = "C:/Users/dpelt/OneDrive - 서울시립대학교/Desktop/Mayson/UOS_graduate/julia"
data = open(@sprintf "%s/diabetes.txt" path)
# data = open("C:/Users/dpelt/OneDrive - 서울시립대학교/Desktop/Mayson/UOS_graduate/julia/diabetes_scaled.txt")
lines = readlines(data)

# data setting
n = length(lines)
p = length(split(lines[1]))
data = zeros(n, p)
for i in 1:length(lines)
    line = split(lines[i])
    data[i, p] = parse(Float64, line[1])
    for j in 2:p
        data[i, (j-1)] = parse(Float64, line[j][3:end])
    end
end
X = data[:, 1:p-1]
y = data[:, end]
p -= 1

# data scaling (?)
for i in 1:p
    X[:, i] = normalize!(X[:, i])
end
norm(X[:, 1])

# parameters
C = 10
tol = 1e-6
max_passes = 1e+2

# initalize
alpha = zeros(n, 1)
# alpha = rand(n, 1)
# b = rand(1)
b = 0
passes = 0

Random.seed!(520)
# main
while passes < max_passes
    num_changed_alphas = 0

    for i in 1:n
        f_i = sum((alpha .* y) .* (X * X[i, :])) + b
        E_i = f_i - y[i]

        if((y[i] * E_i < -tol) & (alpha[i] < C) | (y[i] * E_i > tol) & (alpha[i] > 0))
            j = rand(1:n)
            while true
                j = rand(1:n)
                if j != i
                    break
                end
            end

            f_j = sum((alpha .* y) .* (X * X[j, :])) + b
            E_j = f_j - y[j]

            # save old alpha
            alpha_old_i = alpha[i]
            alpha_old_j = alpha[j]

            if y[i] != y[j]
                L = max(0, alpha_old_j - alpha_old_i)
                H = min(C, C + alpha_old_j - alpha_old_i)
            else 
                L = max(0, alpha_old_j + alpha_old_i - C)
                H = min(C, alpha_old_j + alpha_old_i)
            end
          
            if L == H
                continue
            end

            eta = 2 * dot(X[i, :], X[j, :]) - dot(X[i, :], X[i, :]) - dot(X[j, :], X[j, :])
            if eta >= 0
                continue
            end

            alpha[j] = alpha_old_j - y[i] * (E_i - E_j) / eta

            if abs(alpha[j] - alpha_old_j) < 1e-5
                continue
            end

            alpha[i] = alpha_old_i + y[i] * y[j] * (alpha_old_j - alpha[j])

            b1 = b - E_i - y[i] * (alpha[i] - alpha_old_i) * dot(X[i, :], X[i, :]) - 
                y[j] * (alpha[j] - alpha_old_j) * dot(X[i, :], X[j, :])
    
            b2 = b - E_j - y[i] * (alpha[i] - alpha_old_i) * dot(X[i, :], X[j, :]) - 
                y[j] * (alpha[j] - alpha_old_j) * dot(X[j, :], X[j, :])

            if (0 < alpha[i]) & (alpha[i] < C)
                global b = b1
            elseif (0 < alpha[j]) & (alpha[j] < C)
                global b = b2
            else 
                global b = (b1 + b2) / 2
            end

            num_changed_alphas += 1
        end
    end

    if num_changed_alphas == 0
        global passes += 1
    else
        global passes = 0
    end
end

# result check
result = zeros(n, 1)
for k in 1:n
    result[k, :] = sum((alpha .* y) .* (X * X[k, :])) + b
end
sum(result .>= 0) == n # all predictions are 1