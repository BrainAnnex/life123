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
# ## Reversible Synthesis elementary reaction `A + B <-> C`
# #### taken to equilibrium.
# #### Comparison of 2 approximate solutions, and the exact solutions 

# %% [markdown]
# ### TAGS :  "uniform compartment", "numerical"

# %%
LAST_REVISED = "Oct. 17, 2026"
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
# # 1. Run a simulation with low accuracy

# %% [markdown]
# ### Initialize the System
# Specify the chemicals, the reactions, and the initial concentrations

# %%
# Instantiate the simulator and specify the chemicals
# Here we use the "fast" preset for the variable steps, which leads to fewer steps, but generally less-accurate results
uc_fast = ReactionSimulator(method="forward_euler")    # We call this "uc" because it's a "uniform compartment

# %%
# Elementary reaction A + B <-> C
uc_fast.add_reaction(reactants=["A" , "B"], products="C", 
                     kinetic_parameters={"kF": 5., "kR": 2.})

uc_fast.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
uc_fast.set_conc({"A": 10., "B": 50., "C": 20.})

# %%

# %% [markdown]
# ### Run the reaction

# %%
uc_fast.single_compartment_react(initial_step=0.004, duration=0.02, variable_steps=False) # NOT using variable steps!

# %%
uc_fast.get_history()

# %% [markdown]
# ### Plots changes of concentrations with time

# %%
uc_fast.plot_history(show_intervals=True)

# %%

# %% [markdown]
# # 2. Repeat the simulation with the "Heun" method, which is more sophisticated
# ### and, generally, more accurate

# %%
# Instantiate the simulator and specify the chemicals
uc_slow = ReactionSimulator(method="heun", reaction_registry=uc_fast.get_reactions())
# Re-use the species and reactions of part 1, but now with the "heun" ODE solver

# %%
# Set the initial concentrations of all the chemicals
uc_slow.set_conc({"A": 10., "B": 50., "C": 20.})

# %%
uc_slow.single_compartment_react(initial_step=0.004, duration=0.02, variable_steps=False) # NOT using variable steps!

# %%
uc_slow.get_history()

# %%
uc_slow.plot_history(show_intervals=True)

# %%

# %% [markdown]
# # 3. And, finally, get the EXACT analytic solution

# %%
# Instantiate the simulator and specify the chemicals
uc_exact = ReactionSimulator(reaction_registry=uc_fast.get_reactions(),
                             method="analytic")
# Re-use the species and reactions of part 1, but now with an exact ODE solver

# %%
# Set the initial concentrations of all the chemicals
uc_exact.set_conc({"A": 10., "B": 50., "C": 20.})

# %%
uc_exact.single_compartment_react(initial_step=0.001, duration=0.02, variable_steps=False) 
# Here we use smaller steps to more accurately depict the exact solution in the whole interval

# %%
#uc_exact.get_history() 

# %%
uc_exact.plot_history(show_intervals=True)

# %%

# %% [markdown]
# # 4. Let's compare plots for the concentration of `C`, as a function of time, from the earlier 3 simulations

# %%
p1 = uc_fast.plot_history(species="C", colors=['#F5B914'], title="Forward Euler (simpler)", show=False)
p2 = uc_slow.plot_history(species="C", colors=['#CBE504'], title="Heun (more complex)", show=False)
p3 = uc_exact.plot_history(species="C", colors=['forestgreen'], title="exact", show=False)

# %%
PlotlyHelper.combine_plots([p1, p2, p3], title="Comparison of simulation accuracies")

# %% [markdown]
# ## The comparison reveals that the `Heun` method greatly outperformed the `Forward Euler` one
# in this scenario with large fixed simulation steps 

# %%
# TODO:  rename uc_slow and uc_fast   ; show dotted lines at the steps  ; maybe rename C to P (and also in react_3_b)
