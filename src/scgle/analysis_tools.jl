# src/scgle/analysis_tools.jl

"""
    inv_mobility(t_values::AbstractVector{Float64}, Δζ_values::AbstractVector{Float64})

Calculates the (normalized) inverse mobility by integrating the time-dependent friction Δζ(t).
The integral is performed using the trapezoidal rule.
The result is `1 + ∫Δζ(t)dt`.

# Arguments:
* `t_values::AbstractVector{Float64}`: Vector of time points.
* `Δζ_values::AbstractVector{Float64}`: Vector of Δζ(t) values corresponding to `t_values`.

# Returns:
* `Float64`: The calculated inverse mobility.
"""
function inv_mobility(t, Δζ)
    dt = diff(t)
    db⁻¹ = 0.5*dt.*(Δζ[1:end-1]+Δζ[2:end])
    return 1 + sum(db⁻¹)
end


"""
    get_Δr²(solution_object_main_scgle, solver_options, 
            k_array_full::Vector{Float64}, σ_vector::Vector{Float64}, 
            particle_species_index::Integer, αc_val::Float64;
            MSD0::Float64=0.0, dMSD0::Float64=0.0) # dMSD0 usualmente 6D_short_time

Calcula el Desplazamiento Cuadrático Medio (MSD) para una especie de partícula dada.
Primero extrae la fricción Δζ(t) de `solution_object_main_scgle`, luego construye
un `MSDSCGLEKernel` y resuelve la ecuación de memoria para el MSD.

# Arguments:
* `solution_object_main_scgle`: Solución de la ecuación SCGLE principal (de `BinarySCGLEKernel`).
* `solver_options`: Opciones para el solver de `ModeCouplingTheory.jl` (e.g., `MCTSolverOptions`).
* `k_array_full::Vector{Float64}`: Vector de números de onda (se usa `k_array_full[1]` como k₀).
* `σ_vector::Vector{Float64}`: Vector de diámetros de partícula [σ₁, σ₂].
* `particle_species_index::Integer`: Índice de la especie (1 o 2).
* `αc_val::Float64`: Constante αc para `get_Δζ`.
* `MSD0::Float64=0.0`: Condición inicial para MSD(t=0).
* `dMSD0::Float64=0.0`: Condición inicial para dMSD/dt(t=0). (Para F = MSD, esta es MSD'(0). Si la ecuación es para dMSD/dt, esta es (dMSD/dt)'(0) ).

# Returns:
* `Tuple`: (`t_values`, `MSD_values_species`) - Tiempos y valores del MSD para la especie especificada.
"""
function get_Δr²(sol, solver, k, σ)
    MSD0 = 0.0; dMSD0 = 0.0; α = 0.0; β = 1.0; γ = 0.0; δ = -1.0;
    msdkernel1 = MSDSCGLEKernel(sol, k, σ, 1)
    msdequation1 = MemoryEquation(α, β, γ, δ, MSD0, dMSD0, msdkernel1)
    msdsol1 = solve(msdequation1, solver)
    msdkernel2 = MSDSCGLEKernel(sol, k, σ, 2)
    msdequation2 = MemoryEquation(α, β, γ, δ, MSD0, dMSD0, msdkernel2)
    msdsol2 = solve(msdequation2, solver)
    return msdsol1.t, msdsol1.F, msdsol2.F
end

#########
#   η   #
#########
"""
ΔG = (kBT/60π^2)∫dkk⁴[∂C⋅F]²
"""
function get_ΔG(sol, k_array, Sₖ)
    ΔG₁ = []
    ΔG₂ = []
    t = sol.t
    Nk = div(length(k_array),2)
    Δk = k_array[2] - k_array[1]
    k⁴ = k_array[Nk+1:end-1].^4
    S⁻¹ = inv.(Sₖ[Nk+1:end])
    Cₖ = similar(Sₖ[Nk+1:end])
    δαβ = Matrix(I(2))
    for i in 1:Nk
        Cₖ[i] = (δαβ - S⁻¹[i]) # Eq. (71) Naegele-Banchio
    end
    ∂C = diff(Cₖ)
    F = sol.F
    for it in 1:length(t)
        dG = zeros(2,2)
        for ik in 1:(Nk-1)
            A = ∂C[ik]*F[it][Nk+ik]
            dG += k⁴[ik]*(A*A')
        end
        append!(ΔG₁, dG[1,1]*Δk/(60*π^2))
        append!(ΔG₂, dG[2,2]*Δk/(60*π^2))
    end
    return ΔG₁, ΔG₂
end

function get_τₐ(τ, F)
    idx = sum(F .> exp(-1))
    return τ[idx]
end
