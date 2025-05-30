module ThermodynamicInstabilities

# Es posible que necesites importar funciones específicas de AnaliticalStructureFactors.jl
# si S_RPA_mixture_SquareWell no se importa a nivel del paquete NESCGLE2.jl.
# Ejemplo:
# using AnaliticalStructureFactors: S_RPA_mixture_SquareWell
# O si NESCGLE2.jl ya hace `using AnaliticalStructureFactors`, entonces
# `AnaliticalStructureFactors.S_RPA_mixture_SquareWell` debería estar disponible.

export calculate_M_elements, calculate_fluctuation_angle, 
       calculate_inverse_S_min_eigenvalue, check_gibbs_criterion

"""
    calculate_M_elements(σ_vector::Vector{Float64}, ϕ_vector::Vector{Float64}, 
                         Amplitude_matrix::Matrix{Float64}, λ_matrix::Matrix{Float64}, 
                         composition_vector::Vector{Float64}; k_value::Float64 = 0.0)

Calcula los elementos `M_ρρ`, `M_ρc`, `M_cc` de la matriz de estabilidad termodinámica
basada en la representación de Bhatia-Thornton, a partir de los factores de estructura
parciales en el límite de `k_value` (generalmente `k=0`).

Estos elementos están relacionados con las segundas derivadas de la energía libre de Gibbs.
La matriz `M` es:
`[[M_ρρ, M_ρc], [M_ρc, M_cc]]`

# Argumentos
- `σ_vector::Vector{Float64}`: Vector de diámetros de las partículas de las especies.
- `ϕ_vector::Vector{Float64}`: Vector de fracciones de volumen de las especies.
- `Amplitude_matrix::Matrix{Float64}`: Matriz de amplitudes de la interacción de pozo cuadrado.
- `λ_matrix::Matrix{Float64}`: Matriz de rangos de la interacción de pozo cuadrado.
- `composition_vector::Vector{Float64}`: Vector con las fracciones molares `[x₁, x₂]` de las dos especies.
- `k_value::Float64 = 0.0`: Valor del número de onda para el cual se calculan los factores de estructura. Por defecto es 0.0 para el límite termodinámico.

# Retorna
- `Tuple{Float64, Float64, Float64}`: Una tupla conteniendo `(M_ρρ, M_ρc, M_cc)`.

# Dependencias
- Requiere una función `S_RPA_mixture_SquareWell` (presumiblemente de `AnaliticalStructureFactors.jl`)
  que calcule la matriz de factores de estructura parciales `s_ij(k)`.

# Referencias
- Las fórmulas para `M_ρρ`, `M_cc`, `M_ρc` se basan en la transformación de los
  factores de estructura de Bhatia-Thornton (ver, por ejemplo, las ecuaciones
  2.18, 2.20, 2.21 mencionadas por el usuario).
"""
function calculate_M_elements(s_matrix::Matrix{Float64}, composition_vector::Vector)
    if length(composition_vector) != 2
        error("composition_vector debe tener dos elementos [x₁, x₂].")
    end
    x1 = composition_vector[1]
    x2 = composition_vector[2]

    
    A = inv(s_matrix) # Inversa de la matriz de factores de estructura S_ij(k=0)

    # M_ρρ = (x₁ * A₁₁ + x₂ * A₂₂ + 2 * √(x₁x₂) * A₁₂)
    M_ρρ = (x1 * A[1,1]) + (x2 * A[2,2]) + (2 * sqrt(x1 * x2) * A[1,2])
    
    # M_cc = (x₂ * A₁₁ + x₁ * A₂₂ - 2 * √(x₁x₂) * A₁₂)
    M_cc = (x2 * A[1,1]) + (x1 * A[2,2]) - (2 * sqrt(x1 * x2) * A[1,2])
    
    # M_ρc = √(x₁x₂) * (A₁₁ - A₂₂) - (x₁ - x₂) * A₁₂
    M_ρc = sqrt(x1*x2)*(A[1,1] - A[2,2]) - (x1-x2)*A[1,2]
    
    return M_ρρ, M_ρc, M_cc
end

"""
    calculate_fluctuation_angle(M_ρρ::Float64, M_ρc::Float64, M_cc::Float64)

Calcula el ángulo `α` que caracteriza la dirección de la fluctuación dominante (eigenvector
correspondiente al eigenvalor más pequeño `λ₁`) de la matriz de estabilidad termodinámica `M`.
La matriz `M` es `[[M_ρρ, M_ρc], [M_ρc, M_cc]]`.

El ángulo se define tal que `tan(α) = x₁ρ / x₁c`, donde `x = [x₁ρ, x₁c]` es el eigenvector.
La relación se deriva de `(M_ρρ - λ₁)x₁ρ + M_ρc x₁c = 0`.

# Argumentos
- `M_ρρ::Float64`: Elemento (1,1) de la matriz de estabilidad.
- `M_ρc::Float64`: Elemento (1,2) y (2,1) de la matriz de estabilidad.
- `M_cc::Float64`: Elemento (2,2) de la matriz de estabilidad.

# Retorna
- `Float64`: El ángulo `α` en radianes.

# Notas
- Maneja casos especiales donde el denominador `(M_ρρ - λ₁)` es cercano a cero.
"""
function calculate_fluctuation_angle(M_ρρ::Float64, M_ρc::Float64, M_cc::Float64)
    # λ₁ es el eigenvalor más pequeño de la matriz de estabilidad M
    # λ₁ = 0.5 * (Tr(M) - sqrt(Tr(M)² - 4*det(M)))
    # Tr(M) = M_ρρ + M_cc
    # det(M) = M_ρρ*M_cc - M_ρc²
    # sqrt((M_ρρ + M_cc)² - 4*(M_ρρ*M_cc - M_ρc²)) = sqrt(M_ρρ² + M_cc² + 2M_ρρM_cc - 4M_ρρM_cc + 4M_ρc²)
    # = sqrt(M_ρρ² + M_cc² - 2M_ρρM_cc + 4M_ρc²) = sqrt((M_ρρ - M_cc)² + 4M_ρc²)
    
    discriminant = (M_ρρ - M_cc)^2 + 4 * M_ρc^2
    if discriminant < 0.0 # Debería ser no negativo para M real y simétrica
        # Esto podría ocurrir por errores numéricos si es muy cercano a cero.
        # Asumimos que es cero en tal caso.
        discriminant = 0.0
    end
    
    λ1 = 0.5 * (M_ρρ + M_cc - sqrt(discriminant))
    
    denominator = M_ρρ - λ1
    
    # Precisión de máquina para Float64
    machine_epsilon = eps(Float64)

    if abs(denominator) < machine_epsilon 
        # Caso: M_ρρ - λ₁ ≈ 0
        # Esto implica que M_ρρ es (o está muy cerca de) λ₁.
        # La ecuación de eigenvector es M_ρc * x₁c ≈ 0.
        if abs(M_ρc) < machine_epsilon
            # Caso: M_ρc ≈ 0 también.
            # M es diagonal o casi diagonal, y M_ρρ es el eigenvalor más pequeño.
            # Eigenvector es [1, 0] (o un múltiplo). x₁ρ/x₁c → Inf.
            # Ángulo es π/2 (fluctuación a lo largo de ρ).
            # Si M_ρρ = M_cc y M_ρc = 0, λ1 = M_ρρ. Denominador es 0.
            # Cualquier vector es un eigenvector. Por convención, podemos elegir [1,0] o [0,1].
            # Si M_ρρ < M_cc y M_ρc = 0, λ1 = M_ρρ. Denominador es 0. Eigenvector [1,0]. α = π/2.
            return π/2.0 
        else
            # Caso: M_ρc ≠ 0.
            # Entonces x₁c debe ser 0. Eigenvector es [x₁ρ, 0], ej. [1, 0].
            # tan(α) = x₁ρ/0 → ±Inf.
            # El signo de atan(±Inf) es ±π/2.
            # tanα = -M_ρc / (M_ρρ - λ1)
            # Si M_ρρ - λ1 -> 0⁺, tanα -> -sign(M_ρc) * Inf
            # Si M_ρρ - λ1 -> 0⁻, tanα -> +sign(M_ρc) * Inf
            # atan(-M_ρc/denominator) manejará esto.
            # Sin embargo, para ser explícito: si x₁c = 0, la fluctuación es puramente en ρ.
            # α = ±π/2. El signo depende de la convención de M_ρc.
            # -M_ρc / (un número muy pequeño)
            return atan(-M_ρc, denominator) # atan(y,x) es más robusto para cuadrantes
        end
    else 
        # Caso: M_ρρ - λ₁ ≠ 0 (caso general)
        tanα = -M_ρc / denominator
        return atan(tanα)
    end
end

"""
    calculate_inverse_S_min_eigenvalue(S_matrix::Matrix{Float64})

Calcula el eigenvalor más pequeño (`λ⁻`) de la inversa de una matriz `S` 2x2 dada.
Se asume que `S` es simétrica.

# Argumentos
- `S_matrix::Matrix{Float64}`: Una matriz 2x2, usualmente la matriz de factores de estructura `S_ij(k)`.

# Retorna
- `Float64`: El eigenvalor más pequeño de `inv(S_matrix)`.

# Errores
- Lanza un error si `S_matrix` no es 2x2.
"""
function calculate_inverse_S_min_eigenvalue(S_matrix::Matrix{Float64})
    if size(S_matrix) != (2,2)
        error("La matriz de entrada S_matrix debe ser 2x2 para este cálculo específico de eigenvalor.")
    end

    # Verificar si S_matrix es singular podría ser útil antes de inv()
    # det_S = S_matrix[1,1]*S_matrix[2,2] - S_matrix[1,2]*S_matrix[2,1]
    # if abs(det_S) < eps(Float64)
    #     @warn "La matriz S_matrix es singular o casi singular. La inversa puede ser inestable."
    #     return Inf # O algún otro valor indicador
    # end

    A = inv(S_matrix) # Matriz inversa A = S⁻¹
    
    # El eigenvalor más pequeño de una matriz simétrica A 2x2 es:
    # λ⁻ = 0.5 * ( (A₁₁ + A₂₂) - √((A₁₁ - A₂₂)² + 4*A₁₂²) )
    # Asumiendo A₁₂ = A₂₁
    discriminant = (A[1,1] - A[2,2])^2 + 4*A[1,2]^2 # A[1,2] * A[2,1] si no es simétrica, pero inv(simétrica) es simétrica.
    if discriminant < 0.0
        discriminant = 0.0 # Por errores numéricos
    end

    λ_minus = 0.5 * ( (A[1,1] + A[2,2]) - sqrt(discriminant) )
    
    return λ_minus
end

"""
    check_gibbs_criterion(M_ρρ::Float64, M_ρc::Float64, M_cc::Float64)

Calcula una métrica para la inestabilidad de Gibbs, dada por `M_cc - M_ρc² / M_ρρ`.
Esta cantidad es proporcional al determinante de la matriz de estabilidad `M` (det(M)/M_ρρ).
Un valor negativo (asumiendo `M_ρρ > 0`) indica inestabilidad.

# Argumentos
- `M_ρρ::Float64`: Elemento (1,1) de la matriz de estabilidad.
- `M_ρc::Float64`: Elemento (1,2) y (2,1) de la matriz de estabilidad.
- `M_cc::Float64`: Elemento (2,2) de la matriz de estabilidad.

# Retorna
- `Float64`: El valor del criterio de Gibbs. Puede ser `-Inf` si `M_ρρ` es cero y `M_ρc` no lo es.

# Notas
- Maneja el caso donde `M_ρρ` es cercano a cero para evitar división por cero.
"""
function check_gibbs_criterion(M_ρρ::Float64, M_ρc::Float64, M_cc::Float64)
    machine_epsilon = eps(Float64)

    if abs(M_ρρ) < machine_epsilon
        # Si M_ρρ ≈ 0:
        if abs(M_ρc) < machine_epsilon
            # Si M_ρc ≈ 0 también, entonces det(M) ≈ 0 * M_cc - 0² = 0.
            # El criterio se vuelve indeterminado (0/0) o simplemente M_cc si se interpreta
            # como el límite de la estabilidad de la concentración pura.
            # Si M_ρρ es exactamente 0, la estabilidad depende de M_cc (si M_ρc también es 0).
            return M_cc 
        else
            # Si M_ρc ≠ 0, entonces det(M) ≈ -M_ρc².
            # El criterio M_cc - M_ρc²/M_ρρ tendería a -∞ (si M_ρρ > 0) o +∞ (si M_ρρ < 0),
            # indicando fuerte inestabilidad si M_ρρ se acerca a 0 desde el lado positivo.
            # Un valor grande y negativo es un indicador de inestabilidad.
            return -Inf # Indica inestabilidad fuerte
        end
    end
    return M_cc - (M_ρc^2 / M_ρρ)
end

end # module ThermodynamicInstabilities
