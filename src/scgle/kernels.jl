# src/scgle/kernels.jl

# ===================================
# Core SCGLE Kernel (BinarySCGLEKernel)
# ===================================

struct BinarySCGLEKernel{T1, T2} <: ModeCouplingTheory.MemoryKernel
    k::T1 # Vector of wavenumbers at which the structure factor is known
    λ::T2 # Empirical kernel-weight factor
    cn::T2 # Partial direct correlation function
    S⁻¹nh::T2 # Complement partial direct correlation function
    prefactor::T1 # Escaling factor kernel
end

"""
    BinarySCGLEKernel(ρ, σ, k_array, S_array, S⁻¹_array)
Constructor for the memory function of the form
    K(k, t) = λ(k)Δζ(t)
where
    λₐᵦ(k) = δₐᵦ/(1+(k/kcₐ)²)
and 
    Δζ(t)ₐᵦ = δₐᵦΔζ(t)ᵦ
with
    Δζ(t)ₐ = D⁰ₐ/6π² ∫dk k⁴[Fₛ(t)]ₐₐ[c⋅√n⋅F(t)⋅S⁻¹⋅√n⋅h]ₐₐ

# Arguments:
* `ρ`: vector of number densities for each species
* `σ`: vector of diameters for each species
* `k_array`: vector of wavenumbers at which the structure factor is known
* `S_array`: a `Vector` of Nk `SMatrix`s containing the structure factor of each component at each wave number
* `S⁻¹_array`: a `Vector` of Nk inverse of `SMatrix`s containing the structure factor of each component at each wave number

# Returns:

An instance `k` of `ModeCouplingKernel <: MemoryKernel`, which can be called both in-place and out-of-place:
`kernel = BinarySCGLEKernel(ρ, σ, k, S, S⁻¹)`
`evaluate_kernel!(out, kernel, F, t)`
`out = evaluate_kernel(kernel, F, t)`
"""
function BinarySCGLEKernel(ρ, σ, k_array, S_array, S⁻¹_array)
    Nk = length(k_array)
    Nk2 = div(Nk, 2)
    Δk = k_array[2] - k_array[1]
    #prefactor = Δk*D/(6*π*π)
    prefactor = Δk*(σ.^(-1))/(6*π*π)
    # lambda
    λ = [zeros(2,2) for i in 1:Nk]
    kc = (αc_scgle*2π)*(σ.^(-1))
    # Weights
    U = reshape([1.0, 0.0, 0.0, 1.0], (2,2))
    N = reshape([√(ρ[1]), 0.0, 0.0, √(ρ[2])], (2,2))
    N⁻¹ = inv(N)
    cn = [zeros(2,2) for i in 1:Nk2]
    S⁻¹nh = [zeros(2,2) for i in 1:Nk2]
    for (i, k) in enumerate(k_array[1:Nk2])
        λ[i][1,1] = 1/(1+(k/kc[1])^2)
        λ[i][2,2] = 1/(1+(k/kc[2])^2)
        λ[i+Nk2][1,1] = 1/(1+(k/kc[1])^2)
        λ[i+Nk2][2,2] = 1/(1+(k/kc[2])^2)
        S = S_array[Nk2+i]
        S⁻¹ = S⁻¹_array[Nk2+i]
        cn[i] = N⁻¹*(U - S⁻¹)
        S⁻¹nh[i] = (U - S⁻¹)* N⁻¹
    end
    return BinarySCGLEKernel(k_array, λ, cn, S⁻¹nh, prefactor)
end


function ModeCouplingTheory.evaluate_kernel!(out::Diagonal, 
                                            kernel::BinarySCGLEKernel, 
                                            F, 
                                            t)
    k_array = kernel.k
    Nk = length(k_array)
    Nk2 = div(Nk,2)
    # lambda
    λ = kernel.λ
    # weights
    cn = kernel.cn
    S⁻¹nh = kernel.S⁻¹nh
    # dζ(q) integral
    dζ = [0.0, 0.0]
    for i in 1:Nk2
        k⁴ = k_array[i]^4
        M = cn[i]*F[Nk2+i]*S⁻¹nh[i]
        dζ[1] += F[i][1,1]*M[1,1]*k⁴
        dζ[2] += F[i][2,2]*M[2,2]*k⁴
    end
    dζ[1] *= kernel.prefactor[1]
    dζ[2] *= kernel.prefactor[2]
    for i in 1:Nk
        out.diag[i] = reshape([λ[i][1,1]*dζ[1], 0.0, 0.0, λ[i][2,2]*dζ[2]], (2,2))
    end
end

function ModeCouplingTheory.evaluate_kernel(kernel::BinarySCGLEKernel, F::Vector, t)
    out = Diagonal(similar(F)) # we need it to produce a diagonal matrix
    evaluate_kernel!(out, kernel, F, t) # call the inplace version
    return out
end

"""
    get_Δζ(sol::MemoryEquationSolution, k₀::Float64, σ::Vector{Float64})
Obtains Δζ factor from the the kernel `K` from a `MemoryEquationSolution` object
    K(k, t) = λ(k)Δζ(t)
where
    λₐᵦ(k) = δₐᵦ/(1+(k/kcₐ)²)

Then, 
Δζₐₐ(t) = K(k₀, t)/λₐₐ(kₒ)
# Arguments:
* `sol`: MemoryEquationSolution object
* `k₀`: initial value of the vector of wavenumbers
* `σ`: vector of diameters for each species

# Returns:

Two arrays of the Memory function Δζ
"""
function get_Δζ(sol, k₀::Float64, σ::Vector{Float64})
    K = get_K(sol)
    Δζ₁₁ = []
    Δζ₂₂ = []
    λ₁₁ = 1/((1+(k₀/(αc_scgle*2*π*σ[1]))^2))
    λ₂₂ = 1/((1+(k₀/(αc_scgle*2*π*σ[2]))^2))
    for idx in 1:size(K)[1]
        append!(Δζ₁₁, K[idx][1,1][1,1]/λ₁₁)
        append!(Δζ₂₂, K[idx][1,1][2,2]/λ₂₂)
    end
    return Δζ₁₁, Δζ₂₂
end

# ===================================
# MSD Kernel (MSDSCGLEKernel)
# ===================================

struct MSDSCGLEKernel{TDICT,T} <: MemoryKernel
    tDict::TDICT
    Δζ::T
end

"""
MSDSCGLEKernel(ϕ, D⁰, k_array, Sₖ, sol)

Constructor of a MSDSCGLEKernel. It implements the kernel

K(k,t) = D⁰ / (36πϕ) ∫dq q^4 c(q)^2 F(q,t) Fs(q,t)

where the integration runs from 0 to infinity. F and Fs are the coherent
and incoherent intermediate scattering functions, and must be passed in
as solutions of the corresponding equations.

# Arguments:

* `ϕ`: volume fraction
* `D⁰`: Short times diffusion coefficient
* `k_array`: vector of wavenumbers at which the structure factor is known
* `Sₖ`: vector with the elements of the structure factor 
* `sol`: a solution object of an equation with a SCGLEKernel.

# Returns:

an instance `k` of `MSDSCGLEKernel <: MemoryKernel`, which can be evaluated like:
`k = evaluate_kernel(kernel, F, t)`
"""
function MSDSCGLEKernel(sol, k_array, σ, particle::Integer)
    Δζ₁, Δζ₂ = get_Δζ(sol, k_array[1], σ)
    t = sol.t
    tDict = Dict(zip(t, eachindex(t)))
    kernel = particle == 1 ? MSDSCGLEKernel(tDict, Δζ₁) : MSDSCGLEKernel(tDict, Δζ₂)
    return kernel
end

function evaluate_kernel(kernel::MSDSCGLEKernel, MSD, t)
    it = kernel.tDict[t]
    return kernel.Δζ[it]
end
