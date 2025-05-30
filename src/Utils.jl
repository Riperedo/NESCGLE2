# src/Utils.jl
module Utils

export bisection

"""
    bisection(condition::Function, A::Real, T_target::Real, tolerance::Real; flag = false)

Function to perform a sucessive bisection between two interest points.
It finds the value `Achilles` starting from `A` towards `T_target`
such that `condition(Achilles)` is true, and `condition(Achilles + smallest_step)`
would be false (if searching towards increasing values and condition is for "being below a threshold")
or `condition(Achilles - smallest_step)` would be false (if searching towards decreasing values
and condition is for "being above a threshold").

Essentially, it finds the boundary where the condition changes.
The interpretation of 'A' and 'T_target' (e.g., high_temp, low_temp) depends on how
the `condition` function is defined and the search direction.

# Arguments
- `condition::Function`: Criteria to perform a step. Should return `true` to continue moving `Achilles`.
- `A::Real`: Initial point for `Achilles`.
- `T_target::Real`: The other boundary of the search interval for `Tortoise`.
- `tolerance::Real`: The minimum step size to consider for convergence.

# Keywords
- `flag::Bool = false`: If true, prints internal computing steps.
"""
function bisection(condition::Function, A::Real, T_target::Real, tolerance::Real; 
    flag = false)
    Achilles = A
    Tortoise = T_target 
    δ = 0.5 # Initial step factor
    paso(t, a, Delta) = Delta*(t - a) 

    initial_step_direction = sign(Tortoise - Achilles)
    if initial_step_direction == 0
        return Achilles 
    end

    while true 
        current_step_size = paso(Tortoise, Achilles, δ)

        if abs(current_step_size) < tolerance
            break
        end

        next_Achilles = Achilles + current_step_size

        if condition(next_Achilles) 
            Achilles = next_Achilles
        else
            Tortoise = next_Achilles 
            δ = 0.5 
        end

        if Achilles == Tortoise || sign(Tortoise - Achilles) != initial_step_direction
            break
        end
        if flag println("Achilles: ", Achilles, ", Step: ", current_step_size, ", Next_Achilles: ", next_Achilles, ", Tortoise: ", Tortoise) end
    end
    return Achilles
end

end # module Utils