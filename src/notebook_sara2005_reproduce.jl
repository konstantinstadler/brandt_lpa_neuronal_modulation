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

# ╔═╡ f85c6a11-3c29-431b-9016-7c02c8bc533a
# We first load libraries necessary for this excercise
begin
    import Pkg
	Pkg.activate(mktempdir())

    Pkg.add("DifferentialEquations")
    Pkg.add("Plots")
    Pkg.add("PlutoUI")
	using DifferentialEquations
	using Plots
	using PlutoUI
	plotly()
end

# ╔═╡ 2e8ba40e-a759-11eb-10c8-eb85ce6d59a4
md"# Reproducing the vesicle recycling model from  [Sara et al 2005](https://doi.org/10.1016/j.neuron.2004.12.056)"

# ╔═╡ be5a053c-a2fb-456e-ba4c-830996bdb9a8
md"This notebook implements the differential equations for the model of spontaneous synaptic vesicle recycling. The model is based on single pool with four states and in the article fitted against experimental observation with tagged/dyed vesicles.
"


# ╔═╡ a58c1ee9-8274-4030-b2e2-247592f5d862
md"The original description can be found at [Sara et al 2005. An Isolated Pool of Vesicles Recycles at Rest and Drives Spontaneous Neurotransmission. Neuron 45/4](https://doi.org/10.1016/j.neuron.2004.12.056): Figure 5 and Experimental Procedures (\"Modeling Spontaneous Synaptic Vesicle Recycling\")."

# ╔═╡ 491066b3-600f-4a79-87d0-3e49ff648e41
md"
The model as in the article describes a single pool with four states:

- s0, dye-loaded vesicles
- s1, mobilized vesicles
- s2, empty vesicles
- s3, recycled and mixed empty vesicles

with parameters

- α rate of mobilization 
- β rate of recycling
- δ rate of dye loss; 

Note, that Julia is one (1) indexed, the equations in the article are 0 (zero) indexed:"

# ╔═╡ d129cc09-c4cd-471a-8a8c-0d092b105de3

function sara_spon_recycle!(ds, s, p, tspan)
# Note: Julia is 1 indexed, the formulas in the original publication are 0 indexed
    α, β, δ = p.α, p.β, p.δ 
    ds[1] = -α * s[1] * (s[1] / (s[1] + s[4]))
    ds[2] = -δ * s[2] + α * (s[1] / (s[1] + s[4])) * s[1]
    ds[3] = -β * s[3] + δ*s[2]
    ds[4] = β * s[3]
end


# ╔═╡ 6b242b10-bdec-4bb0-b681-6fd4ea9c62e7
md"Next we define the parameters and initial values for the model. Here we use given values in the model description at the Experimental procedures and Fig 5. 

**Note** that the parameter values (p below) have contradicting values in the article text:

α (rate of mobilization) is set to 0.0008 $s^{-1}$. The motivation given is that this is 1 vesicle release per 120 sec. This, however, would lead to an α of 1/120 = 0.0083 (Caption text of Figure 5 in the article). Also the updated release rate (for 8mM of extracellular Ca) - 1 vesicle every 60 sec (0.0166) - does not correspond to the given value of 0.0016


"

# ╔═╡ 49cf7f94-d8b9-4878-b483-ecffc5f5307f
begin

md"set α (rate of mobilization) to $(@bind alpha_enter NumberField(0.0001:0.0001:0.1,default=0.0008))

set β (rate of recycling) to $(@bind beta_enter NumberField(0.0:0.01:1.0,default=0.5))

set δ (rate of dye loss) to $(@bind sig_enter NumberField(0.0:0.01:2.0,default=1.67))
"

	
end

# ╔═╡ 6541cc0b-3990-450a-84d9-1df196e81e74
begin
	
	
	# paper states a pool of 15 vesicles, but present graphs normalized to 1 - we set all in the dye-loaded state s1 (s0 in the paper)
	s_init = [1.0,0.0,0.0,0.0];
	
	# origianl graph was for up to 1200 sec
	t_span = (0.0,1200.0)
	
	# values from the sliders
	para = (α=alpha_enter, β=beta_enter, δ=sig_enter);
end

# ╔═╡ fef63cff-662d-45dc-975e-486cb954595c
md"Next we setup the ODE solver and solve with the given parameters"

# ╔═╡ 88409be8-bb1a-40a4-9f92-f62132eece7b
prob = ODEProblem(sara_spon_recycle!,s_init,t_span,para);

# ╔═╡ 497d59a4-3e20-4935-91c8-e7d496d0d2f4
sol = solve(prob, Tsit5(), abstol = 1e-9, reltol = 1e-9);

# ╔═╡ ae55cf1f-45bd-4eaa-9a04-6b204d786d8d
begin
	Plots.theme(:juno)
	plot(sol,vars=(0,1), label="s1(t)", margin=5Plots.mm)
	ylabel!("Fraction of loaded pool")
	title!("Loaded pool (reprod. Fig 5D Sara et al 2005)", fontsize=10)
end

# ╔═╡ 10a090cb-4e0d-4070-a6a1-a31da9068d0f
begin
	p2 = plot(sol,vars=(0,2), label="s2(t)", title="mobilized vesicles", xlabel="");
	p3 = plot(sol,vars=(0,3), label="s3(t)", title="empty vesicles", xlabel="");
	p4 = plot(sol,vars=(0,4), label="s4(t)", title="recycled and mixed ves.");
	plot(p2, p3, p4, layout=grid(3,1), margin=5Plots.mm)
end

# ╔═╡ Cell order:
# ╟─2e8ba40e-a759-11eb-10c8-eb85ce6d59a4
# ╟─be5a053c-a2fb-456e-ba4c-830996bdb9a8
# ╟─a58c1ee9-8274-4030-b2e2-247592f5d862
# ╟─f85c6a11-3c29-431b-9016-7c02c8bc533a
# ╟─491066b3-600f-4a79-87d0-3e49ff648e41
# ╟─d129cc09-c4cd-471a-8a8c-0d092b105de3
# ╟─6b242b10-bdec-4bb0-b681-6fd4ea9c62e7
# ╟─49cf7f94-d8b9-4878-b483-ecffc5f5307f
# ╟─6541cc0b-3990-450a-84d9-1df196e81e74
# ╟─fef63cff-662d-45dc-975e-486cb954595c
# ╟─88409be8-bb1a-40a4-9f92-f62132eece7b
# ╟─497d59a4-3e20-4935-91c8-e7d496d0d2f4
# ╟─ae55cf1f-45bd-4eaa-9a04-6b204d786d8d
# ╟─10a090cb-4e0d-4070-a6a1-a31da9068d0f
