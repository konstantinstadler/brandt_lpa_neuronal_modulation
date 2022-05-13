# LPA Vesicle model without endocytosis reduction
# This mirrors the conditions in vivo slices.
# 
# To run in the julia interpreter (after activating the environment - see readme):
# include("./src/lpa_vesicle_model_only_mobilization_increase.jl")

using DifferentialEquations
using DataFrames
using CSV
using Plots

# Theming of plots
# -----------------
Plots.theme(:juno)

# Storing results
# ---------------

const RESULT_DIR = normpath(joinpath(@__FILE__,"..","..", "results"))
mkpath(RESULT_DIR)
# FIG_NAME = "lpa_vesicle_model.png"
# TABLE_NAME = "lpa_vesicle_model.csv"

FIG_NAME = "lpa_vesicle_model_only_mobilization_increase.png"
TABLE_NAME = "lpa_vesicle_model_only_mobilization_increase.csv"

# Parameters
# -----------

# Initial parameters
p = (α=1/120, β=0.5, σ=1.67)
t_exp_duration = (0.0,1200.0)  # in sec
u0 = [1.0,0.0, 0.0]            # particle state at the start

t_lpa_effect = 200
p_lpa_effect = (α=1.5/120, β=0.5, σ=1.67)


function vesicle_model!(du, u, p, t_exp_duration)
    α, β, σ = p.α, p.β, p.σ
    du[1] = -α * u[1] + β * u[3]      # vesicles ready for release
    du[2] = +α * u[1] - σ * u[2]      # vesicles currently merged with the membrane
    du[3] = +σ * u[2] - β * u[3]      # currently being recycled vesicles
end


function first_lpa_time_affect(y, t, integrator)
    return t < t_lpa_effect
end

function first_para_affect!(integrator)
    integrator.p = p_lpa_effect
end


callbacks = CallbackSet(
    ContinuousCallback(first_lpa_time_affect, first_para_affect!))


prob = ODEProblem(ODEFunction(vesicle_model!, syms=[:u1, :u2, :u3]) , u0, t_exp_duration, p)
sol = solve(prob, Tsit5(), abstol = 1e-9, reltol = 1e-9, callback=callbacks)

df_result = DataFrame(sol)
CSV.write(joinpath(RESULT_DIR, TABLE_NAME), df_result, delim='\t')

plot(sol, 
     vars = (0,2),
     xlim = (10,1200),      # allow some time for reaching a steady state
     xlabel = "time [s]",
     ylabel = "Fraction of currently merged vesicles")
savefig(joinpath(RESULT_DIR, FIG_NAME))
# gui()  # to show the plot uncomment



