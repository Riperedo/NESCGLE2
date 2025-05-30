using NESCGLE2
using StaticArrays 
using AnalyticalStructureFactors


ϕ_a = 0.3 # fraccion de volumen de particulas a
ϕ_b = 0.281 # fraccion de volumen de particulas b

delta = 0.8 # razon de diametros de particulas b con respecto a particulas a 

ϕ_vector = [ϕ_a, ϕ_b]
σ_vector = [1.0, delta]  # diametros

# Wave vector variables
Nk = 200
kmax = 15*π; dk = kmax/Nk
k_array = dk*(collect(1:Nk) .- 0.5)

# convert the data to the Vector of SMatrix format
S_array = [@SMatrix(zeros(2, 2)) for i = 1:Nk]
for i = 1:Nk
    S_array[i] = S_HS_VW_mixture(σ_vector, ϕ_vector, k_array[i])
end

# SCGLE dynamics
τ, Δζ₁, Δζ₂, W1, W2, ΔG₁, ΔG₂ = NESCGLE2.SelfConsistentScheme.dynamics(ϕ_vector, σ_vector, k_array, S_array)

save_data("dynamics.dat", [τ Δζ₁ Δζ₂ W1 W2 ΔG₁ ΔG₂], header_lines = ["1 τ\t|2 Δζ1\t|3 Δζ2\t|4 W1\t|5 W2"] )

println()
