# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Reversible Synthesis elementary reaction `A + B <-> P`
# #### taken to equilibrium.
# #### Comparison of 2 numerical solvers ("Heun" and "Forward Euler") against the exact solution.
# We'll use large, fixed-sized, steps for the numerical solvers.   
#
# SEE ALSO: experiment `react_3_b`

# %% [markdown]
# ### TAGS :  "uniform compartment", "numerical"

# %%
LAST_REVISED = "Oct. 18, 2026"
LIFE123_VERSION = "1.0.0rc9"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import pandas as pd
from life123 import check_version, ReactionSimulator,PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# # 1. Run a simulation using the **"Forward Euler"** method (a simple, straighforward method)

# %% [markdown]
# ### Initialize the System
# Specify the chemicals, the reactions, and the initial concentrations

# %%
# Instantiate the simulator
sim_euler = ReactionSimulator(method="forward_euler")

# %%
# Elementary reaction A + B <-> P
sim_euler.add_reaction(reactants=["A" , "B"], products="P", 
                       kinetic_parameters={"kF": 5., "kR": 2.})

sim_euler.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
sim_euler.set_conc({"A": 10., "B": 50., "P": 20.})

# %%

# %% [markdown]
# ### Run the reaction

# %%
sim_euler.single_compartment_react(initial_step=0.004, duration=0.02, variable_steps=False) # NOT using variable steps!

# %%
sim_euler.get_history()

# %% [markdown]
# ### Plots changes of concentrations with time

# %%
sim_euler.plot_history(show_intervals=True)

# %% [markdown]
# The simulator clearly got to ambitious in that 1st step!  (In fact, if you try with a slightly larger step size, it'll give an error about excessive step size)

# %%

# %%

# %% [markdown]
# # 2. Repeat the simulation with the **"Heun"** method, which is more sophisticated
# ### and, generally, more accurate.  
# SAME FIXED STEP SIZE

# %%
# Instantiate a new simulator
sim_heun = ReactionSimulator(method="heun", reaction_registry=sim_euler.get_reactions())
# Re-use the species and reactions of part 1, but now with the "heun" ODE solver

# %%
# Set the initial concentrations of all the chemicals
sim_heun.set_conc({"A": 10., "B": 50., "P": 20.})

# %%
# NOT using variable steps!   Same fixed step as before
sim_heun.single_compartment_react(initial_step=0.004, duration=0.02, variable_steps=False) 

# %%
sim_heun.get_history()

# %%
sim_heun.plot_history(show_intervals=True)

# %%

# %%

# %% [markdown]
# # 3. And, finally, get the EXACT analytic solution

# %%
# Instantiate a new simulator
sim_exact = ReactionSimulator(reaction_registry=sim_euler.get_reactions(),
                              method="analytic")
# Re-use the species and reactions of part 1, but now with an exact ODE solver

# %%
# Set the initial concentrations of all the chemicals
sim_exact.set_conc({"A": 10., "B": 50., "P": 20.})

# %%
sim_exact.single_compartment_react(initial_step=0.0005, duration=0.02, variable_steps=False) 
# Here we use smaller steps to more accurately depict the exact solution in the whole interval

# %%
#sim_exact.get_history() 

# %%
sim_exact.plot_history(show_intervals=True)

# %%

# %%

# %% [markdown]
# # 4. Let's compare plots for the concentration of the product `P`, as a function of time, from the earlier 3 simulations

# %%
p1 = sim_euler.plot_history(species="P", colors=['#F5B914'], title="Forward Euler (simpler)", show_intervals=True, show=False)
p2 = sim_heun.plot_history(species="P",  colors=['#CBE504'], title="Heun (more complex)", show=False)
p3 = sim_exact.plot_history(species="P", colors=['forestgreen'], title="EXACT", show=False)

# %%
PlotlyHelper.combine_plots([p1, p2, p3], title="Comparison of simulation accuracies")

# %% [markdown]
# ## The comparison reveals that the `Heun` method greatly outperformed the `Forward Euler` one
# in this scenario with large, fixed simulation steps.  
# Notice how `Forward Euler` starts out with exactly the correct initial slope (that's exactly how it works) - but that ends up being its downfall when continuing along that direction for an extended time step!

# %%
