### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ ad798920-094d-42d0-9740-76e90ca54a84
begin
    import Pkg
    Pkg.add("DifferentialEquations")
    Pkg.add("Plots")
    Pkg.add("PlutoUI")

	using DifferentialEquations
	using Plots
	using PlutoUI
	plotly()
end

# ╔═╡ cb14eb66-a829-11eb-010a-43bf9cbe798a
md"# Simple vesicle recycling model"

# ╔═╡ 096f1405-280a-41f4-b67b-ce052aba3c24
md"A simple two pool model of spontaneous vesicle recycling. It is a simplified version of the model published in [Sare et al 2005](https://doi.org/10.1016/j.neuron.2004.12.056)"

# ╔═╡ 7b94a7cc-bd31-4129-bed0-2df4e4443bdb
md" We assume two pools:

- u1: Vesicles currently in the resting state
- u2: Vesicles currently activated/merged with the the pre-synaptic membrane

in addition we have two parameters:

- α: activation/exocytosis rate of vesicles in u1
- β: recycling rate from the membrane back to the resting pool"

# ╔═╡ 57e6296c-a95b-437c-9254-9267ef576b90
# Differential equations describing the system
function vesicle_recycle!(du, u, p, t_span)
    α, β = p.α, p.β
    du[1] = -α * u[1] + β * u[2]      # vesicles in the resting pool 
    du[2] = +α * u[1] - β * u[2]      # fused vesicle currently at the membrane
end



# ╔═╡ 74d113d9-e358-482c-bd56-55cc422c13d0
md"Setting the initial state of the system"

# ╔═╡ e1ed9123-715e-4c78-9a50-cde91d7825b8
# We start with all vesicles in the resting pool
u0 = [1.0,0.0]

# ╔═╡ 327d620b-fbad-48c5-9784-e19254c65ea8
# running the model for 1200 secs
t_span = (0.0,1200.0)

# ╔═╡ 6bcc1592-6e2c-41a3-b19e-d59a95ef5245
p = (α=0.008, β=0.5)

# ╔═╡ 820f98f4-9f2b-4139-b184-deb579f15d4e
md"α after LPA application: $(@bind alpha_new NumberField(0.0001:0.0001:0.1, default=0.016))
	
β after LPA application: $(@bind beta_new NumberField(0.0:0.001:1.0, default=0.05))
	
time of LPA application: $(@bind time_new Slider(0:10:1200, default=200, show_value=true))
	"


# ╔═╡ a90ffa7f-4269-464f-a087-3c93a78718ae


# ╔═╡ f5007ad3-2f95-4229-a5ae-8367c03e80c7
begin
	
	function lpa_application(y, t, integrator)
	    return t < time_new
	end
	
	function affect!(integrator)
	    integrator.p = (α=alpha_new, β=beta_new)
	end
	
	callback = ContinuousCallback(lpa_application, affect!)
end

# ╔═╡ f4b08f7a-3857-4792-a32d-4da2f007b1f1
begin
	prob = ODEProblem(vesicle_recycle!,u0,t_span,p);
	sol = solve(prob, Tsit5(), abstol = 1e-9, reltol = 1e-9, callback=callback);
end;

# ╔═╡ 5efc8b40-84f4-446c-a936-a5dac0699ad0
begin
	Plots.theme(:juno)
	plot(sol, margin=5Plots.mm)
	
end

# ╔═╡ Cell order:
# ╟─cb14eb66-a829-11eb-010a-43bf9cbe798a
# ╟─096f1405-280a-41f4-b67b-ce052aba3c24
# ╟─7b94a7cc-bd31-4129-bed0-2df4e4443bdb
# ╟─ad798920-094d-42d0-9740-76e90ca54a84
# ╟─57e6296c-a95b-437c-9254-9267ef576b90
# ╟─74d113d9-e358-482c-bd56-55cc422c13d0
# ╟─e1ed9123-715e-4c78-9a50-cde91d7825b8
# ╟─327d620b-fbad-48c5-9784-e19254c65ea8
# ╟─6bcc1592-6e2c-41a3-b19e-d59a95ef5245
# ╟─f4b08f7a-3857-4792-a32d-4da2f007b1f1
# ╟─820f98f4-9f2b-4139-b184-deb579f15d4e
# ╠═a90ffa7f-4269-464f-a087-3c93a78718ae
# ╟─f5007ad3-2f95-4229-a5ae-8367c03e80c7
# ╟─5efc8b40-84f4-446c-a936-a5dac0699ad0
