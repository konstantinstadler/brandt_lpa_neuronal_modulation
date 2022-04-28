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

# ╔═╡ 9e472ad5-8de4-4b82-b84d-ce01797fcdbd
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

# ╔═╡ e0dc4b68-a8b7-11eb-0e1c-c7a342d92771
md"# Vesicle recycling model with membrane remain factor"

# ╔═╡ 6fc0f5f6-8d71-4a64-bf2c-9bddeee4b41f
md"A simple two pool model of spontaneous vesicle recycling. It is a simplified version of the model published in [Sara et al 2005](https://doi.org/10.1016/j.neuron.2004.12.056).

Here we change the model from Sara 2005 to remove the pool of empty-recycled vesicles.
Instead, recycled vesicles go back into the loaded vesicle pool
"

# ╔═╡ be6da4c0-b437-40c6-9674-ef51d8a66340
md" We assume three pools:

- u1: Vesicles currently in the resting state
- u2: Vesicles currently activated/merged with the the pre-synaptic membrane
- u3: Vesicles currently being recycled after endocytosis

in addition we have three parameters:

- α: activation/exocytosis rate of vesicles in u1
- β: recycling rate from the membrane back to the resting pool
- σ: vesicle endo-cytosis from the memebrane rate

"



# ╔═╡ e9cd9aed-da67-4083-ada6-2773f5ee4e30
function vesicle_recycle!(du, u, p, t_span)
    α, β, σ = p.α, p.β, p.σ
    du[1] = -α * u[1] + β * u[3]      # vesicles in the resting pool 
    du[2] = +α * u[1] - σ * u[2]      # vesicles currently merged with the membrane
    du[3] = +σ * u[2] - β * u[3]      # currently being recycled vesicles
end

# ╔═╡ c5fcdd8e-1042-417e-9ddb-ba10adc6a67b
md"Setting the initial state of the system"

# ╔═╡ 4c5b59c9-c639-40ef-8c14-d61ebbc71262
# We start with all vesicles in the resting pool
u0 = [0.5, 0.0, 0.5]

# ╔═╡ 8e584a9a-8af4-47fc-8f4c-025b526d26d7
# running the model for 1200 secs
t_span = (0.0,1200.0)

# ╔═╡ 14d8b628-6734-4f89-a3ff-79103fe72038
# initial parameters 
p = (α=0.008, β=0.5, σ=1.67)

# ╔═╡ 008c7875-3c93-4dfb-ad19-57b82fff85a2
md"α after LPA application: $(@bind alpha_new NumberField(0.0001:0.0001:0.01, default=0.018))
	
β after LPA application: $(@bind beta_new NumberField(0.0:0.01:1.0, default=0.01))

σ after LPA application: $(@bind sigma_new NumberField(0.0:0.01:3.0, default=1.67))
	
time of LPA application: $(@bind time_new Slider(0:10:1200, default=200, show_value=true))
	"


# ╔═╡ 58273622-33aa-4e99-b8b3-eeca0985c1d9
begin
	
	function lpa_application(y, t, integrator)
	    return t < time_new
	end
	
	function affect!(integrator)
	    integrator.p = (α=alpha_new, β=beta_new, σ=sigma_new)
	end
	
	callback = ContinuousCallback(lpa_application, affect!)
	
	prob = ODEProblem(vesicle_recycle!,u0,t_span,p);
	sol = solve(prob, Tsit5(), abstol = 1e-6, reltol = 1e-6, callback=callback);
	
	Plots.theme(:juno)
	plot_memb = plot(sol, vars=(0,2), 
		title = "Vesicle in the membrane", 
		xlabel="", ylabel="Ves. Fraction", 
		lc="Dark Orchid")
	plot_all = plot(sol, vars=[(0,1),(0,3)], 
		title="Other pools", 
		xlabel="t (sec)", ylabel="Ves. Fraction")
	plot(plot_memb, plot_all, layout=(2,1), margin=5Plots.mm)

end

# ╔═╡ e611cc3b-b831-433e-848b-dfe05fc7faa4


# ╔═╡ 19db778c-629e-4848-8b9b-9f4f4a03a0a2


# ╔═╡ Cell order:
# ╟─e0dc4b68-a8b7-11eb-0e1c-c7a342d92771
# ╟─6fc0f5f6-8d71-4a64-bf2c-9bddeee4b41f
# ╟─be6da4c0-b437-40c6-9674-ef51d8a66340
# ╠═9e472ad5-8de4-4b82-b84d-ce01797fcdbd
# ╠═e9cd9aed-da67-4083-ada6-2773f5ee4e30
# ╟─c5fcdd8e-1042-417e-9ddb-ba10adc6a67b
# ╠═4c5b59c9-c639-40ef-8c14-d61ebbc71262
# ╟─8e584a9a-8af4-47fc-8f4c-025b526d26d7
# ╟─14d8b628-6734-4f89-a3ff-79103fe72038
# ╟─008c7875-3c93-4dfb-ad19-57b82fff85a2
# ╠═58273622-33aa-4e99-b8b3-eeca0985c1d9
# ╠═e611cc3b-b831-433e-848b-dfe05fc7faa4
# ╠═19db778c-629e-4848-8b9b-9f4f4a03a0a2
