function path_init(
    X, y;
    Aeq = zeros(0, size(X, 2)),
    beq = zeros(size(Aeq, 1)),
    solver = ECOSSolver(maxit=10e8, verbose=0)
    ) 

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
    #solver = ECOSSolver(maxit=10e8, verbose=0)
    solve!(problem, solver)
    target_value = target.valuemeho
    rho_max = target_value[1, :]
    lambda_max = target_value[(p+2):end, :]

    # violation check
    idx1 = findall(abs.((-X' * y + rho_max .* ones(p)) - Aeq' * lambda_max) .<= 1e-4)
    idx2 = findall(abs.((-X' * y - rho_max .* ones(p)) - Aeq' * lambda_max) .<= 1e-4)
    active_set = union(idx1, idx2)

    # Aeq rank 확인하는 부분 추가하면 좋을듯!

    return rho_max, lambda_max, active_set
end

#path_init(X, y, Aeq, beq)