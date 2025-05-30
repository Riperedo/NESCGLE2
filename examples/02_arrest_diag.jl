using NESCGLE2
using StaticArrays 
using AnalyticalStructureFactors


#################
# Initial setup #
#################

const ϕmax = 0.7
const ϕmin = 1e-6

# space arrays
phi_b = collect(0.01:0.01:0.15)

# saving arrays
phi_a_save = Float64[]
phi_b_save = Float64[]

# diameters
const σ_aa = 1.0
const σ_bb = 0.2

# objects to call the main functions
const σ = [σ_aa, σ_bb]

# wave vector
Nk = 200; kmax = 15*π; dk = kmax/Nk
k_array = dk*(collect(1:Nk) .- 0.5)

U = reshape([1.0, 0.0, 0.0, 1.0], (2,2))
sₖ = [convert(SMatrix{2, 2, Float64}, U) for i = 1:Nk]

# main loop
function buscar_estado(estado::String)
	for idx in eachindex(phi_b)
	    ϕb = phi_b[idx]
        printstyled("Looking for $(estado)\n", color=:magenta)
        # primero limitamos la busqueda a ϕmax
        """
        Evaluamos Gammas, salvo que se encuentre en una inestabilida termodinamica.
        Buscamos que se evalue si `sistema == estado`
        """
        function condition(ϕa::Float64)
            # initial condition
            sistema = "Dump"

            if ϕa + ϕb > ϕmax
                sistema = "Glass"
            else
                # si la fraccion de volumen total es menor que ϕmax continuamos
                ϕ_vector = [ϕa, ϕb]
                # Calcular fracciones molares x_i = ρ_i / ρ_total
                ρ_vector = phi_to_rho_mixture(ϕ_vector, σ)
                ρ_total_mixture = sum(ρ_vector)
                x_vector = ρ_vector ./ ρ_total_mixture
                x1 = x_vector[1]
                x2 = x_vector[2]

                printstyled("ϕ = $(ϕ_vector)\n", color =:yellow)
                # primero evaluamos que el sistema no este en una region espinodal
                S0 = S_HS_VW_mixture(σ, ϕ_vector, 0.0)
                A = inv(S0)
                M_ρρ = x1 * A[1,1] + x2 * A[2,2] + 2 * sqrt(x1 * x2) * A[1,2]
                if M_ρρ < 0 # (ρkTχ)⁻¹ > 0
                    sistema = "Glass" # no es estrictamente Glass solo es para indicar que el estado es arrestado
                    printstyled("Thermodynamic instability\n", color=:red )
                else
                    # Ahora que sabemos que no estamos en una espinodal evaluamos el factor de estructura estatico
                    # convert the data to the Vector of SMatrix format
                    Sₖ = [@SMatrix(zeros(2, 2)) for i = 1:Nk]
                    for (i, k_wavevector) in enumerate(k_array)
                        Sₖ[i] = S_HS_VW_mixture(σ, ϕ_vector, k_wavevector)
                        # En caso de hacer algun valor negativo en la estructura estamos en presencia de algun 
                        # tipo de inestabilidad, por lo que no necesitamos evaluar mas.
                        if (Sₖ[i][1,1] < 0) || (Sₖ[i][2,2] < 0) 
                            sistema = "Glass" # una vez mas, no es necesariamente Glass, es solo para indicar un estado arrestado
                            printstyled("Thermodynamic instability\n", color=:red )
                            break
                        end # end if
                    end # end for
                    # Ahora que estamos seguros que no hay ningun tipo de inestabilidad procedemos
                    # a evaluar gamma
                    if sistema == "Dump"
                        # update the long arrays
                        k = vcat(k_array, k_array)
                        S = vcat(sₖ, Sₖ)
                        S⁻¹ = inv.(S)                        
                        kernel = BinarySCGLEKernel(ρ_vector, σ, k, S, S⁻¹);
                        printstyled("Looking for $(estado)\n", color=:magenta)
                        #iteraciones, gammas, sistema = Asymptotic(kernel, S, flag = true)
                        _, _, sistema, _ = Asymptotics(kernel, S, flag = true)
                    end # end if
                end # end if
            end # end if
            # devolvemos el estado del sistema
            return sistema == estado
        end # end function
        if estado == "Fluid"
            ϕa = bisection(condition, ϕmin, ϕmax, 1e-5; flag = false)
        else
            ϕa = bisection(condition, ϕmax, ϕmin, 1e-5; flag = false)
        end # end if
        append!(phi_a_save, ϕa)
        append!(phi_b_save, ϕb)

        save_data("Kinetic_state_$(estado)_delta_$(σ_bb/σ_aa).dat", [phi_a_save phi_b_save], 
                    header_lines = ["Columnas: ϕa, ϕb"])
	end
end

buscar_estado("Glass")
#buscar_estado("Fluid")