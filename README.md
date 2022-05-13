# Lysophosphatidic Acid selectively modulates excitatory Transmission in hippocampal Neurons

Neuronal model to investigate the neuromodulatory effects of LPA.


## Model description


For modelling LPA effects on the level of the synapse we adapted the single-pool, four-state model developed by [Sara 2005](https://www.sciencedirect.com/science/article/pii/S0896627305000693?via%3Dihub) which describes the dynamic within a synapses observed with dye-loaded vesicles. 

We excluded the parameter describing the dye loss dynamics, which resulted in a new single-pool three-state model with:

- u<sub>1</sub>: Fraction of vesicles ready for release
- u<sub>2</sub>: Fraction of vesicles currently activated/merged with the pre-synaptic membrane
- u<sub>3</sub>: Fraction of vesicles currently being recycled after endocytosis (empty vesicles not yet in the pool of vesicles that can be released)

The dynamics of the model are described with three parameters:

- α: rate of activation/exocytosis 
- β: recycling rate (after endocytosis back to the pool of vesicles ready for release)
- σ: vesicle endocytosis after release

These were set to the values obtained by Sara et al 2005: α=0.008 s<sup>-1</sup>, β=0.5 s<sup>-1</sup>, σ=1.67 s<sup>-1</sup>, and the dynamics are described with:

- du<sub>1</sub> = -α * u<sub>1</sub> + β * u<sub>3</sub>      
- du<sub>2</sub> = +α * u<sub>1</sub> - σ * u<sub>2</sub>      
- du<sub>3</sub> = +σ * u<sub>2</sub> - β * u<sub>3</sub>      


We implemented this model in the scientific programming language [Julia](https://arxiv.org/abs/1209.5145). 
The full source code script to run the model and generate the data/figures is available at the [authors github repository](https://github.com/konstantinstadler/brandt_lpa_neuronal_modulation) .


## Results

Contrary to expectations the increase in pre-synaptic Ca<sub>i</sub> concentration did not led to an increase in the mEPSC frequency.
We hypothesize that this could be explained by a secondary effect of LPA disturbing vesicle recycling in cultured neurons.
To explain this mechanism, we adopted a single-pool pre-synaptic model ([Sara 2005](https://www.sciencedirect.com/science/article/pii/S0896627305000693?via%3Dihub)), which allows to vary the rate of exocytosis α (activation/release of vesicles), endocytosis σ and the recycling rate β (refilling of vesicles and transport to the pool of vescicles ready for release).

We modelled the expected increase in mEPSC frequency by increasing the exocytosis rate α by a factor of 1.5. 
Subsequently, we matched the observed rate of mEPSC reduction by reducing the recycling rate β by a factor of 0.01.
Nevertheless, these model parameters would still result in a transient increase of merged membrane vesicles (not shown). 
However, if the near collapse of the recycling rate happens before the increase of the exocytosis (here shown with a one minute delay) it would mask the transient increase in vesicle release. Such a modelled dynamic matches the observed LPA effect in cultured neurons.


## Running the model

The model was implemented in [Julia](https://julialang.org/) using the [DifferentialEquations](https://juliapackages.com/p/differentialequations) packages. 

To run the model:

- [download and install Julia](https://julialang.org/downloads/)
- download this repository (pressing the button 'Code' - either clone via git or just download/extract the zip file).
- go to the downloaded and extracted directory and start julia (command line: `julia`)
- in Julia activate the environment with 
    - press `]` (to get into the Julia package manager)
    - type `activate .`
    - type `instantiate`
    - the first time it might take some minutes for installing/pre-compiling all packages, it is faster the second time
    - press backspace to go back to the standard Julia prompt
    - run the script with `include("./src/lpa_vesicle_model.jl")`
- this will produce results in the /results folder
    - lpa_vesicle_model.csv : numerical outcome of the model run
    - lpa_vesicle_model.png : visualization of the fraction of merged particles over time (corresponding to u(2) in the csv file)

To also show the output as graph during the run, uncomment `#gui()` at the last line of the script.

The model was developed and tested with Julia version 1.7.2 in Linux 5.17.5 (Pop!_OS 22.04) on a   
Dell 5820 16 core Intel(R) Xeon(R) W-2245 CPU @ 3.90GHz (128GB RAM). Usual runtime is below 5 sec. 

