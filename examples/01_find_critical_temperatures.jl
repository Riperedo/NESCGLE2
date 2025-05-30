# examples/full_sweep_thermodynamic_instabilities.jl
#
# Este script de ejemplo demuestra cómo usar el paquete NESCGLE2.jl para
# realizar un barrido de parámetros (fracciones de volumen) y calcular
# las temperaturas críticas asociadas a inestabilidades termodinámicas
# en una mezcla coloidal binaria. Los resultados se guardan en un archivo.

using NESCGLE2
using AnalyticalStructureFactors # Para S_RPA_mixture_SquareWell y save_data

# --- 1. Definición de Funciones Auxiliares ---
# (Si calculate_mole_fractions no está en NESCGLE2.Utils, definirla aquí)
# Asumiremos que NESCGLE2.Utils.calculate_mole_fractions existe o la defines aquí:
function calculate_mole_fractions_local(ϕ_vector::Vector{Float64}, σ_vector::Vector{Float64})
    if length(ϕ_vector) != 2 || length(σ_vector) != 2
        error("ϕ_vector y σ_vector deben tener dos elementos.")
    end
    if any(σ_vector .<= 0)
        error("Los diámetros de las partículas (σ_vector) deben ser positivos.")
    end
    particle_volumes = (π/6.0) .* (σ_vector.^3)
    ρ_mixture = ϕ_vector ./ particle_volumes
    ρ_total_mixture = sum(ρ_mixture)
    if ρ_total_mixture == 0.0
        return [0.0, 0.0], 0.0 
    end
    x_vector = ρ_mixture ./ ρ_total_mixture
    return x_vector, ρ_total_mixture
end


println("--- Ejemplo de Cálculo de Inestabilidades Termodinámicas (Barrido Completo) con NESCGLE2.jl ---")

# --- 2. Definición de Parámetros del Sistema ---
println("\n[Definiendo parámetros del sistema]")

# Diámetros de las partículas
const σ_aa = 1.0
const σ_bb = 1.0
const σ_vector_params = [σ_aa, σ_bb] # Renombrado para evitar conflicto en bucle

# Rangos de la interacción de pozo cuadrado
const σ_ab_interaction_range_factor_params = 0.5 * (σ_aa + σ_bb) 
const λ_range_aa_params = 1.5 * σ_aa 
const λ_range_bb_params = 1.5 * σ_bb
const λ_range_ab_params = 1.5 * σ_ab_interaction_range_factor_params

const interaction_range_matrix_params = [λ_range_aa_params  λ_range_ab_params; # Renombrado
                                         λ_range_ab_params  λ_range_bb_params]

# Profundidades de pozo de referencia
const u_aa_ref_params = 1.0 # Renombrado
const u_bb_ref_params = 1.0 # Renombrado

# Límite termodinámico (número de onda k -> 0)
const k_thermodynamic_limit_params = 0.0 # Renombrado

# Rangos para el barrido de fracciones de volumen
phi_a_range = collect(0.01:0.05:0.45) # Reducido el paso para un ejemplo más rápido
phi_b_range = collect(0.01:0.05:0.45) # Puedes ajustarlo a 0.01 si prefieres
println("Rango para ϕ_a: ", phi_a_range)
println("Rango para ϕ_b: ", phi_b_range)

# Arrays para almacenar los resultados del barrido
phi_a_results = Float64[]
phi_b_results = Float64[]
tc_eigenvector_results = Float64[]
alpha_at_tc_results = Float64[]
tc_compressibility_results = Float64[]
tc_gibbs_results = Float64[]
rho_total_results = Float64[]
x1_results = Float64[]

println("\n[Iniciando barrido de parámetros...]")

# --- 3. Bucle Principal de Barrido ---
for ϕ_a_loop in phi_a_range
    for ϕ_b_loop in phi_b_range
        print("\rProcesando: (ϕ_a, ϕ_b) = (", round(ϕ_a_loop, digits=2), ", ", round(ϕ_b_loop, digits=2), 
              ")")

        current_ϕ_vector = [ϕ_a_loop, ϕ_b_loop]
        
        # Calcular fracciones molares para el punto actual
        # Usar la función local o la del paquete si existe:
        # current_x_vector, current_ρ_total = NESCGLE2.Utils.calculate_mole_fractions(current_ϕ_vector, σ_vector_params)
        current_x_vector, current_ρ_total = calculate_mole_fractions_local(current_ϕ_vector, σ_vector_params)

        if current_ρ_total == 0.0 # Omitir si la densidad es cero (e.g. ambas phi son 0)
            # (Aunque el bucle empieza en 0.01, esta es una salvaguarda)
            println("\nDensidad total cero para ϕ=", current_ϕ_vector, ", omitiendo.")
            continue
        end
        
        # --- Definición de Funciones de Criterio (cierres sobre variables del bucle) ---
        function compresibility_instability_crit_loop(T::Float64)
            amplitude_matrix = [T/u_aa_ref_params  Inf; Inf  T/u_bb_ref_params]
            s_matrix = AnalyticalStructureFactors.S_RPA_mixture_SquareWell(
                σ_vector_params, current_ϕ_vector, k_thermodynamic_limit_params,
                amplitude_matrix, interaction_range_matrix_params
            )
            m_ρρ, _, _ = NESCGLE2.ThermodynamicInstabilities.calculate_M_elements(s_matrix, current_x_vector)
            return m_ρρ > 0
        end

        function gibbs_instability_crit_loop(T::Float64)
            amplitude_matrix = [T/u_aa_ref_params  Inf; Inf  T/u_bb_ref_params]
            s_matrix = AnalyticalStructureFactors.S_RPA_mixture_SquareWell(
                σ_vector_params, current_ϕ_vector, k_thermodynamic_limit_params,
                amplitude_matrix, interaction_range_matrix_params
            )
            m_ρρ, m_ρc, m_cc = NESCGLE2.ThermodynamicInstabilities.calculate_M_elements(s_matrix, current_x_vector)
            gibbs_metric = NESCGLE2.ThermodynamicInstabilities.check_gibbs_criterion(m_ρρ, m_ρc, m_cc)
            return gibbs_metric > 0
        end

        function eigenvector_instability_crit_loop(T::Float64)
            amplitude_matrix = [T/u_aa_ref_params  Inf; Inf  T/u_bb_ref_params]
            s_matrix = AnalyticalStructureFactors.S_RPA_mixture_SquareWell(
                σ_vector_params, current_ϕ_vector, k_thermodynamic_limit_params,
                amplitude_matrix, interaction_range_matrix_params
            )
            min_eigenvalue = NESCGLE2.ThermodynamicInstabilities.calculate_inverse_S_min_eigenvalue(s_matrix)
            return min_eigenvalue > 0
        end

        # --- Búsqueda de Temperaturas Críticas ---
        # Usar NESCGLE2.Utils.bisection o NESCGLE2.bisection
        T_c_eigen = NESCGLE2.Utils.bisection(eigenvector_instability_crit_loop, 10.0, 1e-5, 1e-5)
        T_c_comp = NESCGLE2.Utils.bisection(compresibility_instability_crit_loop, 10.0, 1e-5, 1e-5)
        T_c_gib = NESCGLE2.Utils.bisection(gibbs_instability_crit_loop, 10.0, 1e-5, 1e-5)

        # --- Cálculo del Ángulo de Fluctuación en T_c_eigen ---
        amplitude_at_Tc_eigen = [T_c_eigen/u_aa_ref_params  Inf; Inf T_c_eigen/u_bb_ref_params]
        s_matrix_at_Tc = AnalyticalStructureFactors.S_RPA_mixture_SquareWell(
            σ_vector_params, current_ϕ_vector, k_thermodynamic_limit_params,
            amplitude_at_Tc_eigen, interaction_range_matrix_params
        )
        m_ρρ_tc, m_ρc_tc, m_cc_tc = NESCGLE2.ThermodynamicInstabilities.calculate_M_elements(
            s_matrix_at_Tc, current_x_vector
        )
        alpha_val = NESCGLE2.ThermodynamicInstabilities.calculate_fluctuation_angle(m_ρρ_tc, m_ρc_tc, m_cc_tc)

        # --- Almacenar Resultados ---
        push!(phi_a_results, ϕ_a_loop)
        push!(phi_b_results, ϕ_b_loop)
        push!(tc_eigenvector_results, T_c_eigen)
        push!(alpha_at_tc_results, alpha_val)
        push!(tc_compressibility_results, T_c_comp)
        push!(tc_gibbs_results, T_c_gib)
        push!(rho_total_results, current_ρ_total)
        push!(x1_results, current_x_vector[1])
    end
end
println("\nBarrido completado.")

# --- 4. Guardar Resultados ---
println("\n[Guardando resultados...]")
output_filename = "spinodal_data_NESCGLE2_example.dat"
header_info = ["Resultados del cálculo de inestabilidades termodinámicas con NESCGLE2.jl",
               "Parámetros del sistema: σ_aa=$(σ_aa), σ_bb=$(σ_bb), λ_aa_range=$(λ_range_aa_params), etc.",
               "Columnas: phi_a | phi_b | Tc_eigenvector | alpha_at_Tc | Tc_compressibility | Tc_Gibbs | rho_total | x1"]

data_to_save = hcat(phi_a_results, phi_b_results, tc_eigenvector_results, 
                    alpha_at_tc_results, tc_compressibility_results, 
                    tc_gibbs_results, rho_total_results, x1_results)

# Usar save_data de AnalyticalStructureFactors
AnalyticalStructureFactors.save_data(output_filename, data_to_save, header_lines=header_info)
println("Resultados guardados en: ", output_filename)

println("\n--- Fin del Ejemplo (Barrido Completo) ---")

# Para ejecutar este ejemplo:
# 1. Asegúrate de tener NESCGLE2.jl y AnalyticalStructureFactors.jl instalados y accesibles.
#    (especialmente que NESCGLE2.Utils.bisection y AnalyticalStructureFactors.save_data estén disponibles).
# 2. Guarda este script como (por ejemplo) examples/run_full_sweep_example.jl
# 3. Desde la REPL de Julia, navega al directorio de ejemplos y ejecuta:
#    include("run_full_sweep_example.jl")