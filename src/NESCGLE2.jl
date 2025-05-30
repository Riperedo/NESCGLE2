# src/NESCGLE2.jl
module NESCGLE2

# Importar dependencias necesarias a nivel de paquete
using AnalyticalStructureFactors
using ModeCouplingTheory

# Incluir el módulo de utilidades
include("Utils.jl")
using .Utils # Para hacer accesibles las funciones de Utils dentro de NESCGLE2
export bisection # Opcional: re-exportar bisection

# Incluir y exportar funcionalidades de los submódulos
include("ThermodynamicInstabilities.jl")
using .ThermodynamicInstabilities
# export ... funciones importantes de ThermodynamicInstabilities ...

include("SelfConsistentScheme.jl")
using .SelfConsistentScheme
export dynamics
export BinarySCGLEKernel


include("AsymptoticProperties.jl")
using .AsymptoticProperties
export Asymptotics
export get_γ
export get_non_ergodic_params


# También puedes re-exportar selectivamente funciones de SelfConsistentScheme 
# si quieres que estén en el espacio de nombres de NESCGLE2 directamente:
# export BinarySCGLEKernel, get_Δr², etc.

# ...

end # module NESCGLE2