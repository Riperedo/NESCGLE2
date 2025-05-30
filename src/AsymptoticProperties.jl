module AsymptoticProperties

struct BinarySCGLEKernel{T1, T2}
    k::T1 # Vector of wavenumbers at which the structure factor is known
    λ::T2 # Empirical kernel-weight factor
    cn::T2 # Partial direct correlation function
    S⁻¹nh::T2 # Complement partial direct correlation function
    prefactor::T1 # Escaling factor kernel
end


"""
Calcula el determinante de una matriz 2x2.
M = [m11 m12; m21 m22]
"""
@inline function det_2x2(M)
    return M[1,1] * M[2,2] - M[1,2] * M[2,1]
end

"""
Calcula la inversa de una matriz 2x2.
Devuelve la matriz inversa y un booleano indicando si la inversión fue exitosa.
"""
@inline function mat_inv_2x2(M)
    d = det_2x2(M)
    # Usar una tolerancia pequeña para la comparación con cero
    if abs(d) < 1e-100 # Un valor muy pequeño, ajustar si es necesario
        #println("Warning: Determinant near zero: ", d)
        return M, false # Devuelve la original o una matriz de error, y false
    end
    inv_d = 1.0 / d
    return [M[2,2]*inv_d  -M[1,2]*inv_d;
            -M[2,1]*inv_d  M[1,1]*inv_d], true
end

"""
Multiplica dos matrices 2x2: A * B.
A = [a11 a12; a21 a22], B = [b11 b12; b21 b22]
"""
@inline function mat_mul_2x2(A, B)
    return [A[1,1]*B[1,1] + A[1,2]*B[2,1]  A[1,1]*B[1,2] + A[1,2]*B[2,2];
            A[2,1]*B[1,1] + A[2,2]*B[2,1]  A[2,1]*B[1,2] + A[2,2]*B[2,2]]
end

"""
Suma dos matrices 2x2: A + B.
"""
@inline function mat_add_2x2(A, B)
    return [A[1,1]+B[1,1] A[1,2]+B[1,2];
            A[2,1]+B[2,1] A[2,2]+B[2,2]]
end

"""
Multiplica una matriz 2x2 por un escalar: s * M.
"""
@inline function scalar_mul_mat_2x2(s::Float64, M)
    return [s*M[1,1] s*M[1,2];
            s*M[2,1] s*M[2,2]]
end


# --- Lógica principal ---

"""
Calcula el kernel asintótico (RHS de la Eq. C4, que es 1/γ).
γ_current_iterate es la estimación actual de γ = diag(γ₁, γ₂).
"""
function asymptotic_kernel(kernel, S_struct_factor, γ_current_iterate)
    # Asegurar que γ_current_iterate tenga valores máximos para evitar problemas numéricos si divergen mucho.
    # Esta es la γ física, no su inversa.
    g11 = min(γ_current_iterate[1,1], 1e15)
    g22 = min(γ_current_iterate[2,2], 1e15)
    
    # γ_term_in_eq es la γ que aparece en las ecuaciones (λ + k²γ)⁻¹
    # Esta es la γ física.
    γ_term_in_eq = [g11 0.0; 0.0 g22]

    k_values = kernel.k
    # Nk define sobre cuántos k-puntos se integra (según la lógica original)
    Nk = div(length(k_values), 2) 
    
    # Acceso a los campos del kernel
    λ_k_array = kernel.λ
    cn_k_array = kernel.cn
    S_inv_nh_k_array = kernel.S⁻¹nh

    # Acumuladores para Δζ∞₁ y Δζ∞₂ (que son 1/γ₁ y 1/γ₂)
    delta_zeta_inf_species = [0.0, 0.0]

    for i in 1:Nk # Itera sobre la primera mitad de k_values según la lógica original
        k_val = k_values[i]
        k_sq = k_val^2
        k_pow4 = k_sq^2 # k⁴

        # Matriz λ para el k actual
        λ_i = λ_k_array[i]
        # Matriz S(k) para el k actual (indexación original S[Nk+i])
        Sk_i = S_struct_factor[Nk+i] 
        # Matriz cn para el k actual
        cn_i = cn_k_array[i]
        # Matriz S⁻¹nh para el k actual
        S_inv_nh_i = S_inv_nh_k_array[i]

        # Calcular ψs = λᵢ * (λᵢ + k²γ)⁻¹
        term_inv_ψs = mat_add_2x2(λ_i, scalar_mul_mat_2x2(k_sq, γ_term_in_eq))
        M1_inv, success_m1 = mat_inv_2x2(term_inv_ψs)
        if !success_m1
            #println("Advertencia: Falla en inversión para ψs en k=", k_val)
            continue # O manejar el error de otra forma
        end
        ψs = mat_mul_2x2(λ_i, M1_inv)

        # Calcular ψ = Sᵢλᵢ * (Sᵢλᵢ + k²γ)⁻¹
        Sλ = mat_mul_2x2(Sk_i, λ_i)
        term_inv_ψ = mat_add_2x2(Sλ, scalar_mul_mat_2x2(k_sq, γ_term_in_eq))
        M2_inv, success_m2 = mat_inv_2x2(term_inv_ψ)
        if !success_m2
            #println("Advertencia: Falla en inversión para ψ en k=", k_val)
            continue # O manejar el error de otra forma
        end
        ψ = mat_mul_2x2(Sλ, M2_inv)
        
        # Calcular cnψnh = cnᵢ * ψ * Skᵢ * S⁻¹nhᵢ (según la derivación original del usuario)
        # La fórmula original en el comentario del usuario era cn[i]*ψ*S[Nk+i]*S⁻¹nh[i]
        # S[Nk+i] es Sk_i
        temp_prod = mat_mul_2x2(cn_i, ψ)
        temp_prod = mat_mul_2x2(temp_prod, Sk_i)
        cnψnh = mat_mul_2x2(temp_prod, S_inv_nh_i)
        
        # Acumular para Δζ∞ᵢ (que es 1/γᵢ)
        # El k⁴ viene de k² implícito en la fórmula de SCGLE y k² de d³k -> 4πk²dk
        # El integrando es ψs[α,α] * {c√n Sλ(Sλ+k²γ)⁻¹√n h}[α,α]
        # que corresponde a ψs[α,α] * cnψnh[α,α] según la estructura del código.
        delta_zeta_inf_species[1] += k_pow4 * ψs[1,1] * cnψnh[1,1]
        delta_zeta_inf_species[2] += k_pow4 * ψs[2,2] * cnψnh[2,2]
    end

    # Prefactor de la integral: Δk / (6π²)
    # Asume espaciado uniforme para Δk basado en los primeros puntos de k_values
    Δk = (Nk > 1 && length(k_values) >= 2) ? (k_values[2] - k_values[1]) : 1.0
    if Nk == 1 && length(k_values) == 1 # Caso de un solo punto k
         Δk = k_values[1] # O alguna otra convención para un solo punto
    elseif length(k_values) < 2 && Nk > 0 # Caso borde
        Δk = 1.0 # Evitar error, aunque la integral sería cuestionable
    end

    prefactor_val = Δk / (6.0 * π^2)

    delta_zeta_inf_species[1] *= prefactor_val
    delta_zeta_inf_species[2] *= prefactor_val
    
    # Devuelve diag(Δζ∞₁, Δζ∞₂) que es diag(1/γ₁, 1/γ₂)
    return [delta_zeta_inf_species[1] 0.0; 0.0 delta_zeta_inf_species[2]]
end

"""
Clasifica el estado del sistema y determina qué especies están arrestadas.
γ_solution es la γ actual = diag(γ₁, γ₂).
relative_change_matrix es (γ_prev - γ_curr) * inv(γ_curr).
"""
function selector(γ_solution, relative_change_matrix, error_relativo::Float64, infinito::Float64)
    is_converged_finite = [false, false] # Convergió a un valor finito
    is_diverged_infinite = [false, false] # Divergió (γ es grande)

    for α in 1:2
        # Comprueba si γ divergió (se volvió muy grande)
        if γ_solution[α,α] >= infinito
            is_diverged_infinite[α] = true
        # Si no divergió, comprueba si convergió a un valor finito
        elseif abs(relative_change_matrix[α,α]) <= error_relativo
            is_converged_finite[α] = true
        end
        # Si no es ninguno de los anteriores, necesita más iteraciones (sigue siendo false, false)
    end

    s1_arrested = is_converged_finite[1]
    s1_fluid = is_diverged_infinite[1]
    s2_arrested = is_converged_finite[2]
    s2_fluid = is_diverged_infinite[2]

    if s1_fluid && s2_fluid
        return "Fluid", "Both species are fluid"
    elseif s1_arrested && s2_arrested
        return "Glass", "Both species are arrested (Double Glass)"
    elseif s1_arrested && s2_fluid
        return "Mixed", "Species 1 arrested, Species 2 fluid (Single Glass)"
    elseif s1_fluid && s2_arrested
        return "Mixed", "Species 1 fluid, Species 2 arrested (Single Glass)"
    else # Cualquier otro caso (uno o ambos no decididos)
        return "Dump", "Iteration has not converged to a final state"
    end
end


"""
Función principal para el cálculo asintótico.
S_struct_factor es el S(k) de entrada.
α_damp es el factor de amortiguación.
"""
function Asymptotics(kernel, S_struct_factor; 
                    flag::Bool = false, 
                    It_MAX::Int = 2000, 
                    initial_gamma_val::Float64 = 0.00001,
                    error_relativo_conv::Float64 = 0.000001,
                    gamma_infinito_threshold::Float64 = 1e15,
                    α_damp::Float64 = 0.5) # Factor de amortiguación
    
    # Configuración de la impresión de progreso
    decimo = It_MAX > 0 ? div(It_MAX, 50) : 1
    if flag
        printstyled("|-------------------------|------------------------| <- 100%\n", color=:green)
        printstyled("|", color=:green)
    end
    
    # Estado final e información
    system_state = "Dump"
    arrested_info = "None"
    
    # Inicialización de γ = diag(γ₁, γ₂)
    # γ_current es la γ física que se itera
    γ_current = [initial_gamma_val 0.0; 0.0 initial_gamma_val]
    γ_previous = copy(γ_current)

    # Ciclo principal de iteración
    it = 0 # Contador de iteraciones efectivas
    for current_loop_iter in 1:It_MAX
        it = current_loop_iter
        if flag && decimo > 0 && current_loop_iter % decimo == 0
            printstyled("#", color=:green)
        end

        # Calcula Δζ∞ = diag(1/γ₁, 1/γ₂) usando la γ_current
        delta_zeta_inf_matrix = asymptotic_kernel(kernel, S_struct_factor, γ_current)

        # Calcula la nueva γ candidata: γ_candidate = inv(Δζ∞)
        γ_candidate, is_invertible = mat_inv_2x2(delta_zeta_inf_matrix)
        
        if !is_invertible
            # Esto sucede si alguna Δζ∞ es cero o muy pequeña, lo que significa que γ es infinito.
            # Marcar las γ correspondientes como infinitas.
            #println("Advertencia: Matriz Δζ∞ no invertible en la iteración ", it)
            # Si delta_zeta_inf_matrix[1,1] es ~0, γ_candidate[1,1] debería ser infinito.
            # Si delta_zeta_inf_matrix[2,2] es ~0, γ_candidate[2,2] debería ser infinito.
            γ_candidate[1,1] = abs(delta_zeta_inf_matrix[1,1]) < 1e-100 ? gamma_infinito_threshold : (1.0/delta_zeta_inf_matrix[1,1])
            γ_candidate[2,2] = abs(delta_zeta_inf_matrix[2,2]) < 1e-100 ? gamma_infinito_threshold : (1.0/delta_zeta_inf_matrix[2,2])
            γ_candidate[1,2] = 0.0 # Asegurar que sea diagonal
            γ_candidate[2,1] = 0.0
        end
        
        # Aplicar amortiguación
        γ_new_damped_11 = (1.0 - α_damp) * γ_current[1,1] + α_damp * γ_candidate[1,1]
        γ_new_damped_22 = (1.0 - α_damp) * γ_current[2,2] + α_damp * γ_candidate[2,2]
        
        γ_previous = copy(γ_current)
        γ_current = [γ_new_damped_11 0.0; 0.0 γ_new_damped_22]

        # Calcular convergencia (cambio relativo)
        # Evitar división por cero si γ_current es muy pequeño, aunque no debería serlo si initial_gamma_val > 0
        rel_change_1 = abs(γ_current[1,1]) > 1e-20 ? (γ_previous[1,1] - γ_current[1,1]) / γ_current[1,1] : (γ_previous[1,1] - γ_current[1,1])
        rel_change_2 = abs(γ_current[2,2]) > 1e-20 ? (γ_previous[2,2] - γ_current[2,2]) / γ_current[2,2] : (γ_previous[2,2] - γ_current[2,2])
        convergencia_matrix = [rel_change_1 0.0; 0.0 rel_change_2]
        
        # Clasificar el estado
        system_state, arrested_info = selector(γ_current, convergencia_matrix, error_relativo_conv, gamma_infinito_threshold)
        
        if system_state != "Dump"
            if flag 
                printstyled(" -> ", system_state, " (", arrested_info, ")", color=:light_cyan)
            end
            break # Salir del bucle si se alcanzó un estado final
        end
    end # Fin del bucle de iteraciones

    if flag
        # Rellenar la barra de progreso si terminó antes
        if it < It_MAX && decimo > 0
            remaining_hashes = div(It_MAX - it, decimo)
            for _ in 1:(remaining_hashes > 0 ? remaining_hashes-1 : 0) printstyled(" ", color=:green) end # -1 por el # ya impreso o para ajustar
        end
        printstyled(" Done!\n", color=:green)
    end
    
    if it == It_MAX && system_state == "Dump"
        arrested_info = "Max iterations reached, convergence failed."
        if flag
             println("Advertencia: Se alcanzó el máximo de iteraciones sin convergencia clara.")
        end
    end

    return it, γ_current, system_state, arrested_info
end

function get_γ(kernel, S_struct_factor; It_MAX = 1000, α_damp = 0.5)
    _, γ, _, _ = Asymptotic(kernel, S_struct_factor, flag=false, It_MAX=It_MAX, α_damp=α_damp)
    return γ[1,1], γ[2,2]
end


function get_non_ergodic_params(kernel, k_array, S_array; It_MAX = 1000, α_damp = 0.5)
    γ1, γ2 = get_γ(kernel, S_struct_factor, It_MAX=It_MAX, α_damp=α_damp)
    Nk = length(k_array)
    ψ = [zeros((2,2)) for _ in 1:Nk]
    U = [1 0: 0 1]
    for i in 1:Nk
        k = k_array[i]
        λ_i = kernel.λ[i]
        Sk = S_array[i]
        λSk⁻¹ = mat_inv_2x2(mat_mul_2x2(λ_i, Sk))
        k²γ = scalar_mul_mat_2x2(k^2, [γ1 0; 0 γ2])
        ψ[i] = mat_inv_2x2(U + mat_mul_2x2(k²γ, λSk⁻¹))
    end
    return ψ, γ1, γ2
end


export Asymptotics, get_γ, get_non_ergodic_params

end