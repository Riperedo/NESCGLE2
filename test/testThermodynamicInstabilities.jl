using Test
using NESCGLE2 # Assuming this brings ThermodynamicInstabilities and its exports
# If functions from ThermodynamicInstabilities are not exported by NESCGLE2,
# you might need: using NESCGLE2.ThermodynamicInstabilities
# We also need access to S_RPA_mixture_SquareWell from AnalyticalStructureFactors
using AnalyticalStructureFactors # Assuming S_RPA_mixture_SquareWell is exported by it or NESCGLE2 re-exports it.
import AnalyticalStructureFactors.S_RPA_mixture_SquareWell

@testset "Thermodynamic Instability Calculations" begin

    @testset "Bisection Function Tests" begin
        greater_than_5(x) = x > 5.0
        # Bisection finds the largest value x starting from A (0.0) moving towards T_target (10.0)
        # for which condition(x) is true, before crossing T_target.
        # If condition(A) is false, it might not move as expected unless logic is adapted.
        # The provided bisection logic: Achilles = A. If condition(Achilles+step) is true, Achilles += step.
        # So if greater_than_5(0.0) is false, it needs to step until it's true.
        # The current bisection seems to find the boundary from A towards T_target where condition *stops* being true,
        # or rather, the last point where it *was* true.

        # Test case 1: Find boundary where x > 5, searching from 0 to 10.
        # Achilles starts at 0. condition(0) is false. Tortoise is 10.
        # It should find a value close to 5.0 from the right if condition is for stability (e.g. >0)
        # The user's bisection: A=initial_high_T, T_target=initial_low_T.
        # It finds the lowest T (Achilles) for which condition(T) is true.
        @test bisection(x -> x > 5.0, 10.0, 0.0, 1e-6) ≈ 5.0 atol=1e-5 
        
        # Test case 2: Find boundary where x < 3, searching from 0 to 10 (T_target).
        # Condition: x < 3.0. A=0.0, T_target=10.0.
        # Achilles starts at 0. condition(0) is true.
        # It should find a value close to 3.0.
        @test bisection(x -> x < 3.0, 0.0, 10.0, 1e-6) ≈ 3.0 atol=1e-5
    end

    @testset "Instability Criteria for a Specific Point" begin
        # --- Definición de Parámetros para un caso de prueba ---
        σ_aa_test = 1.0
        σ_bb_test = 1.0
        σ_vector_test = [σ_aa_test, σ_bb_test]

        λ_aa_val_test = 1.5 * σ_aa_test
        λ_bb_val_test = 1.5 * σ_bb_test
        λ_ab_val_test = 1.5 * σ_aa_test # Asumiendo σ_ab = σ_aa para este ejemplo
        
        λ_matrix_test = [λ_aa_val_test  λ_ab_val_test;
                               λ_ab_val_test  λ_bb_val_test]

        u_aa_test = 1.0 # Profundidad del pozo para la especie a (unidades de energía)
        u_bb_test = 1.0 # Profundidad del pozo para la especie b (unidades de energía)
        # u_ab se maneja a través de Inf en Amplitude_matrix, según el script original.

        k_wavevector_test = 0.0 # Límite termodinámico

        # Elegir un punto específico (ϕ_a, ϕ_b)
        ϕ_a_val_test = 0.1
        ϕ_b_val_test = 0.1
        ϕ_vector_val_test = [ϕ_a_val_test, ϕ_b_val_test]

        # Calcular x_vector (fracciones molares)
        ρ_mixture_val_test = phi_to_rho_mixture(ϕ_vector_val_test, σ_vector_test)
        ρ_total_mixture_val_test = sum(ρ_mixture_val_test)
        x_vector_val_test = ρ_mixture_val_test ./ ρ_total_mixture_val_test
        
        # --- Definición de funciones de criterio para la bisección ---
        # Estas funciones capturan las constantes y parámetros del test desde el ámbito externo.

        function Compresibility_instability_crit(T_val::Float64)
            # Amplitude_matrix construida según el script original del usuario.
            # Se asume que S_RPA_mixture_SquareWell y calculate_M_elements la interpretan correctamente.
            current_Amplitude_matrix = [T_val/u_aa_test  Inf; 
                                        Inf            T_val/u_bb_test]
            s_matrix = S_RPA_mixture_SquareWell(
                σ_vector_test, ϕ_vector_val_test, k_wavevector_test, 
                current_Amplitude_matrix, λ_matrix_test
            )
            
            m_ρρ, _, _ = NESCGLE2.ThermodynamicInstabilities.calculate_M_elements(s_matrix, x_vector_val_test
            )
            return m_ρρ > 0 # Condición de estabilidad: (ρkTχ)⁻¹ > 0
        end

        function Gibbs_instability_crit(T_val::Float64)
            current_Amplitude_matrix = [T_val/u_aa_test  Inf; 
                                        Inf            T_val/u_bb_test]
            s_matrix = S_RPA_mixture_SquareWell(
                σ_vector_test, ϕ_vector_val_test, k_wavevector_test, 
                current_Amplitude_matrix, λ_matrix_test
            )
            
            m_ρρ, m_ρc, m_cc = NESCGLE2.ThermodynamicInstabilities.calculate_M_elements(s_matrix, x_vector_val_test)
            gibbs_metric = NESCGLE2.ThermodynamicInstabilities.check_gibbs_criterion(m_ρρ, m_ρc, m_cc)
            return gibbs_metric > 0 # Condición de estabilidad
        end

        function eigenvector_criteria_crit(T_val::Float64)
            current_Amplitude_matrix = [T_val/u_aa_test  Inf; 
                                        Inf            T_val/u_bb_test]
            
            # Necesitamos S_RPA_mixture_SquareWell de AnalyticalStructureFactors
            s_ij_matrix = S_RPA_mixture_SquareWell(
                σ_vector_test, ϕ_vector_val_test, k_wavevector_test, 
                current_Amplitude_matrix, λ_matrix_test
            )
            min_eigenvalue_inv_S = NESCGLE2.ThermodynamicInstabilities.calculate_inverse_S_min_eigenvalue(s_ij_matrix)
            return min_eigenvalue_inv_S > 0 # Condición de estabilidad
        end

        # --- Ejecutar bisección y pruebas ---
        # Los rangos de bisección [10.0, 1e-5] buscan desde T alta a T baja.
        # Encuentran la T más baja (Achilles) para la cual la condición de estabilidad es verdadera.

        @testset "Bisection for Eigenvector Instability" begin
            T_eigenvector = bisection(eigenvector_criteria_crit, 10.0, 1e-5, 1e-5)
            @test T_eigenvector isa Float64
            @test 1e-5 <= T_eigenvector <= 10.0
            # Ejemplo: @test T_eigenvector ≈ expected_T_eigenvector atol=1e-4
            
            # Calcular ángulo alfa a esta temperatura
            final_Amplitude_matrix = [T_eigenvector/u_aa_test  Inf; Inf  T_eigenvector/u_bb_test]
            s_matrix = S_RPA_mixture_SquareWell(
                σ_vector_test, ϕ_vector_val_test, k_wavevector_test, 
                final_Amplitude_matrix, λ_matrix_test
            )
            m_ρρ_final, m_ρc_final, m_cc_final = NESCGLE2.ThermodynamicInstabilities.calculate_M_elements(s_matrix, x_vector_val_test
            )
            alpha_val = NESCGLE2.ThermodynamicInstabilities.calculate_fluctuation_angle(m_ρρ_final, m_ρc_final, m_cc_final)
            @test alpha_val isa Float64
            @test -π/2 <= alpha_val <= π/2
        end
        
        @testset "Bisection for Compressibility Instability" begin
            T_compressibility = bisection(Compresibility_instability_crit, 10.0, 1e-5, 1e-5)
            @test T_compressibility isa Float64
            @test 1e-5 <= T_compressibility <= 10.0
            # Ejemplo: @test T_compressibility ≈ expected_T_compressibility atol=1e-4
        end
        
        @testset "Bisection for Gibbs Instability" begin
            T_gibbs = bisection(Gibbs_instability_crit, 10.0, 1e-5, 1e-5)
            @test T_gibbs isa Float64
            @test 1e-5 <= T_gibbs <= 10.0
            # Ejemplo: @test T_gibbs ≈ expected_T_gibbs atol=1e-4
        end
    end
end
