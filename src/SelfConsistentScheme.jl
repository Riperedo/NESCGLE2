# src/SelfConsistentScheme.jl
module SelfConsistentScheme

# --- Dependencias Externas ---
# Estas son dependencias que deben estar en el Project.toml de NESCGLE2.jl
using StaticArrays
using LinearAlgebra # Para Diagonal, inv, I (aunque I puede ser @SMatrix [1.0 0; 0 1.0])
using DelimitedFiles # Lo tenías en tu script original, puede ser para cargar datos estáticos en ejemplos/scripts.

# --- Dependencia del Paquete Principal o Hermano ---
# Asumimos que ModeCouplingTheory.jl es una dependencia directa de NESCGLE2.jl
using ModeCouplingTheory 

# --- Importaciones Específicas para Extender o Usar Tipos ---
# Es buena práctica importar explícitamente lo que se extiende o usa intensivamente.
import .ModeCouplingTheory: MemoryKernel, evaluate_kernel!, evaluate_kernel, get_K, solve, MemoryEquation
# Si usas MCTSolverOptions como un tipo explícito en las funciones:
# import .ModeCouplingTheory: MCTSolverOptions 

# --- Constantes Físicas Compartidas para este Módulo ---
"""
Constante física αc utilizada en la definición de los factores de vértice λ(k)
dentro del esquema SCGLE. Su valor es aproximadamente 1.305.
kc_α = αc_scgle * 2π / σ_α
"""
const αc_scgle = 1.305

# --- Inclusión de los Archivos Internos del Módulo ---
# El orden puede importar si hay interdependencias en la definición de tipos vs funciones.
# kernels.jl define los structs (BinarySCGLEKernel, MSDSCGLEKernel).
# analysis_tools.jl define funciones que pueden operar sobre estos o sus resultados.
# Con el constructor de MSDSCGLEKernel modificado, no hay dependencia circular de get_Δζ.
include("scgle/kernels.jl")
include("scgle/analysis_tools.jl")

# --- API Pública del Módulo SelfConsistentScheme ---
# Exportar los structs de los kernels
export BinarySCGLEKernel
export MSDSCGLEKernel

# Exportar las funciones de análisis y cálculo
export get_Δζ
export inv_mobility
export get_Δr² # Esta función ahora está en analysis_tools.jl
export get_ΔG
export get_τₐ

# Las funciones `evaluate_kernel!` y `evaluate_kernel` son métodos que extienden
# funciones de `ModeCouplingTheory.jl`. No necesitan ser exportadas desde aquí
# para que funcionen con el sistema de despacho múltiple de Julia, siempre y cuando
# el usuario también tenga `using ModeCouplingTheory`.
# Esta función podría ir al final de src/SelfConsistentScheme.jl
# o en un archivo src/scgle/scgle_solver_setup.jl que se incluye en SelfConsistentScheme.jl

function dynamics(ϕ, σ, k_array, S_array)
    ρ = (6/π)*(ϕ./(σ.^3))
    # wavelength array length
    Nk = length(k_array)

    # Usefull unit arrays
    U = reshape([1.0, 0.0, 0.0, 1.0], (2,2))
    sₖ = [convert(SMatrix{2, 2, Float64}, U) for _ = 1:Nk]

    # concatenated arrays
    k = vcat(k_array, k_array)
    S = vcat(sₖ, S_array)
    S⁻¹ = inv.(S)

    # Initial conditions
    F₀ = copy(S)
    ∂F0 = [@SMatrix(zeros(2,2)) for _ in 1:2*Nk]
    α = 0.0
    β = 1.0
    γ = similar(S).*0.0
    D = reshape([1/σ[1], 0, 0, 1/σ[2]], (2,2))
    for i in 1:Nk
        γ[i] = k[i]^2*D*S⁻¹[i]
        γ[Nk+i] = k[i]^2*D*S⁻¹[Nk+i]
    end
    δ = @SMatrix zeros(2,2)

    # main process
    kernel = BinarySCGLEKernel(ρ, σ, k, S, S⁻¹);
    equation = MemoryEquation(α, β, γ, δ, S, ∂F0, kernel);
    solver = TimeDoublingSolver(verbose=false, N=16, tolerance=10^-8, max_iterations=10^8, t_max=10.0^12)
    sol = solve(equation, solver);

    # post process
    τ = ModeCouplingTheory.get_t(sol)
    Δζ₁, Δζ₂ = get_Δζ(sol, k_array[1], σ)
    b⁻¹₁ = inv_mobility(τ, Δζ₁)
    b⁻¹₂ = inv_mobility(τ, Δζ₂)
    println("b⁻¹ = [$b⁻¹₁, $b⁻¹₂]")
    println("u = [0.0, 0.0]")
    ΔG₁, ΔG₂ = get_ΔG(sol, k_array, S_array)

    τ, W1, W2 = get_Δr²(sol, solver, k, σ)

    return τ, Δζ₁, Δζ₂, W1, W2, ΔG₁, ΔG₂
end

export dynamics

end # module SelfConsistentScheme