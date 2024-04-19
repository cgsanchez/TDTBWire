const FERMI_DIRAC_EXPONENT_LIMIT = 100.0

function fermidirac(ϵ, μ, β)
    expn = β * (ϵ - μ)
    if expn < FERMI_DIRAC_EXPONENT_LIMIT
        return 1.0
    elseif expn > FERMI_DIRAC_EXPONENT_LIMIT
        return 0.0
    else
        return 1.0 / (exp(expn) + 1.0)
    end
end

function fermilevel(βb, N, evals)
    f(μ) = N / 2 - sum(fermidirac.(evals, μ, βb))
    return find_zero(f, (minimum(evals), maximum(evals)))
end

function canonicaldm(H, T)
    βb = 1.0 / (BOLTZMAN_K * T)
    N = size(H, 1)
    F = eigen(H)
    μ = fermilevel(βb, N, F.values)
    ff = Diagonal(fermidirac.(F.values, μ, βb))
    return F.vectors * (ff * F.vectors'), μ
end