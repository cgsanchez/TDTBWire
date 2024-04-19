# STCH parameters from
#
# A. P. Sutton , T. N. Todorov , M. J. Cawkwell & J. Hoekstra (2001) 
# A simple model of atomic interactions in noble metals based explicitly on electronic
# structure, Philosophical Magazine A, 81:7, 1833-1848, DOI: 10.1080/01418610108216639

# SILVER
# const ϵ = 0.0025495
# const a = 4.09
# const p = 12
# const q = 3
# const c = 401.18
# # const ν = 0.36292
# const ν = 1.0
# const U = 6.26
# const element = :Ag
# const Z = 1

# COPPER
# const ϵ = 0.012611
# const a = 3.61
# const p = 9
# const q = 3
# const c = 112.35
# const ν = 0.24304
# const U = 0.0
# const Z = 1

# GOLD
const ϵ = 0.0078680
const a = 4.08
const p = 11
const q = 4
const c = 139.07
const ν = 0.36361
const U = 7.0
const Z = 1

cutoff(x) = 1.0 - 1.0 / (1.0 + exp(-CUTTOF_DECAY * x))

#hmn(r) = - (ϵ*c/2) * (a / r) ^ q * cutoff(r-1.12*a)

function hmn(r)
    if r ≤ 1.12a
        return -((ϵ * c) / 2) * (a / r)^q
    else
        return 0.0
    end
end

vrep(r) = ϵ * (a / r)^p * cutoff(r - 1.12 * a)

const dnn = a * sqrt(2) / 2
const β = hmn(dnn)
