# Scratchpad - final models are in teh other jl files (as Pluto notebooks)


using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DifferentialEquations

using Plots
Plots.theme(:juno)

# Here we change the model from Sara 2005 to remove the pool of empty-recycled vesicles.
# Instead, recycled vesicles go back into the loaded vesicle pool

function vesicle_recycle!(du, u, p, t_span)
    α, β, σ = p.α, p.β, p.σ
    du[1] = -α * u[1] + β * u[3]      # vesicles in the resting pool 
    du[2] = +α * u[1] - σ * u[2]      # vesicles currently merged with the membrane
    du[3] = +σ * u[2] - β * u[3]      # currently being recycled vesicles
end

p = (α=1/120, β=0.5, σ=1.67)

delta_time_diff = 0

for delta_time_diff in [0, 10, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400]

    function first_lpa_time_affect(y, t, integrator)
        return t < 200
    end

    function first_para_affect!(integrator)
        integrator.p = (α=1/120, β=0.005, σ=1.67)
    end


    function second_lpa_time_affect(y, t, integrator)
        return t < 200 + delta_time_diff
    end

    function second_para_affect!(integrator)
        integrator.p = (α=1.5/120, β=0.005, σ=1.67)
    end

    callbacks = CallbackSet(
        ContinuousCallback(first_lpa_time_affect, first_para_affect!),
        ContinuousCallback(second_lpa_time_affect, second_para_affect!))


    # paper only states a pool of 15 vesicles, here we put all in the dye loaded state s1 (s0 in the paper)
    u0 = [1.0,0.0, 0.0]

    # original graph was for up to 1200 sec
    t_span = (0.0,1200.0)


    prob = ODEProblem(vesicle_recycle!,u0,t_span,p)
    sol = solve(prob, Tsit5(), abstol = 1e-9, reltol = 1e-9, callback=callbacks)

    filename = "effect_time_diff_" * string(delta_time_diff) * "ms.png"
    plot(sol, vars=(0,2))

    # plot(sol)
    savefig(filename)
end



