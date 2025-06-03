# NESCGLE2

`NESCGLE2.jl` es un paquete de Julia diseñado para el estudio teórico y numérico de **mezclas coloidales binarias (bidispersas)**.
Este paquete implementa herramientas basadas principalmente en la Teoría Autoconsistente Generalizada de Langevin (SCGLE, por sus siglas en inglés: Self-Consistent Generalized Langevin Equation) para analizar diversas propiedades de estos sistemas.

## Características Principales

El paquete permite realizar los siguientes tipos de análisis:

1.  **Inestabilidades Termodinámicas**:
    *   Cálculo de temperaturas críticas asociadas a la descomposición espinodal.
    *   Evaluación de criterios de estabilidad termodinámica (e.g., compresibilidad, criterio de Gibbs, eigenvalores de la matriz de susceptibilidad inversa).
    *   Soporte para potenciales de interacción como el de Pozo Cuadrado (Square-Well) mediante el uso de factores de estructura RPA (Random Phase Approximation).

2.  **Diagramas de Arresto Cinético**:
    *   Determinación de las líneas de transición entre estados fluidos y estados cinéticamente arrestados (vidrio coloidal).
    *   Aplicación de la teoría SCGLE para predecir el arresto dinámico en mezclas, comúnmente para sistemas de Esferas Duras (Hard-Spheres).

3.  **Dinámica de Equilibrio**:
    *   Cálculo de propiedades dinámicas en equilibrio, tales como:
        *   Tiempos de relajación estructural.
        *   Coeficientes de fricción dependientes de la frecuencia.
        *   Otros parámetros dinámicos relevantes derivados del esquema autoconsistente de la SCGLE.

## Ejemplos

La carpeta `examples/` contiene scripts que demuestran el uso del paquete para:
*   `01_find_critical_temperatures.jl`: Realizar un barrido de parámetros para encontrar temperaturas críticas de inestabilidad termodinámica.
*   `02_arrest_diag.jl`: Construir un diagrama de fases cinético (fluido-vidrio) para una mezcla binaria.
*   `03_equilibrium_dynamics.jl`: Calcular la dinámica de equilibrio para un punto de estado específico de una mezcla binaria.

## Referencias
*   Todo