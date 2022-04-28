# Lysophosphatidic Acid selectively modulates excitatory Transmission in hippocampal Neurons

Neuronal model to investigate the neuromodulatory effects of LPA.


## Model description


For modelling LPA effects on the level of the synapse we adapted the single-pool models developed by 

[Sara 2005](https://www.sciencedirect.com/science/article/pii/S0896627305000693?via%3Dihub) 

The original model describes the dynamic within a synapses observed with dye-loaded vesicles and consists of four states. We excluded the parameter describing the dye loss dynamics, which resulted in a three-state model with:

- u1: Fraction of vesicles currently in the resting state
- u2: Fraction of vesicles currently activated/merged with the the pre-synaptic membrane
- u3: Fraction of vesicles currently being recycled after endocytosis (empty vesicles not yet in the resting state)

The dynamics of the model are described with three parameters:

- α: rate of activation/exocytosis 
- β: recycling rate (from the membrane back to resting pool)
- σ: vesicle endocytosis after release
These were set to the values obtained by Sara et al 2005: α=0.008 s-1, β=0.5 s-1, σ=1.67 s-1 and the dynamics are described with:

du[1] = -α * u[1] + β * u[3]      
du[2] = +α * u[1] - σ * u[2]      
du[3] = +σ * u[2] - β * u[3]      


We implemented this model in the scientific programming language [Julia](https://arxiv.org/abs/1209.5145). The full source code and script to generate the figures is available at FILL IN.


## Results

Contrary to expectations the increase in pre-synaptic Ca<sub>i</sub> concentration did not led to an increase in the mEPSC frequency.
We hypothesize that this could be explained by a secondary effect of LPA disturbing vesicle recycling in cultured neurons.
To elucidate this mechanism, we developed a three state neuronal model (adpated from [Sara 2005](https://www.sciencedirect.com/science/article/pii/S0896627305000693?via%3Dihub)), which allows to vary the rate of exocytosis α (activation/release of vesicles in the resting pool), endocytosis σ and the recycling rate β (refilling of vesicles and transport to the resting pool).

We modelled the expected increase in mEPSC frequency by increasing the exocytosis rate α by a factor of 1.5. 
Subsequently, we matched the observed rate of mEPSC reduction by reducing the recycling rate β by a factor of 0.01.
Nevertheless, these model parameters would still result in a transient increase of merged membrane vesicles (not shown/SI). 
However, if the near collapse of the recycling rate happens shortly before the increase of the exocytosis it would mask the transient increase in vesicle release. Such a modelled dynamic matches the observed LPA effect in cultured neurons.


## Running the model

The model was implemented in [Julia](https://julialang.org/) using the [DifferentialEquations](https://juliapackages.com/p/differentialequations) packages. 

To run the model:

- [Download and install Julia](https://julialang.org/downloads/)
- Download this repository (pressing the button 'Code' - either via git or just the zip file.
- run the script from the command line with `julia --project ./src/lpa_vesicle_model.jl`
- this will produce a png figure in the the current directory with the results
- the first time it might take some minutes for installing/pre-compiling all packages, it is faster the second time


If you want to run the script from within the julia interpreter:
- go to the downloaded and extracted directory and start julia (command line)
- activate the environment with 
    - press `]` (to get into the Julia package manager)
    - type `activate .`
    - type `instantiate`
    - press backspace to go back to the standard Julia prompt
    - `include("./src/lpa_vesicle_model.jl")`

When run as is, the script re-produces the two results in the /results folder.

- lpa_vesicle_model.csv : numerical outcome of the model run
- lpa_vesicle_model.png : visualization of the fraction of merged particles over time

To also show the output as graph during the run, uncomment `#gui()` at the last line of the script.

