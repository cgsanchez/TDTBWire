function DOS(ϵ, evals, Γ = 0.01)
    -1.0 * imag(sum(map(x -> 1.0 / (ϵ - x + 1.0im * Γ), evals))) / pi
end

function dodos()
    mol = readxyz("ag.13.xyz")
    H0 = Matrix{Float64}(undef, mol.natoms, mol.natoms)
    d = Matrix{Float64}(undef, mol.natoms, mol.natoms)
    distances!(d, mol)
    H0!(H0, d, mol)
    ρ, fermilevel = canonicaldm(H0, 100.0)
    print("fermilevel ", fermilevel, "\n")
    F = eigen(H0)
    dos = map(ϵ -> DOS(ϵ, F.values), -20:0.01:20)
    plot(-20:0.01:20, dos)
end

