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
LAST_REVISED = "Jan. 4, 2026"
LIFE123_VERSION = "1.0.0rc7"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import pandas as pd
from life123 import check_version, UniformCompartment,PlotlyHelper

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
uc_fast = UniformCompartment(preset="fast")    

# %%
# Elementary reaction A + B <-> C
uc_fast.add_reaction(reactants=["A" , "B"], products="C", kF=5., kR=2.)

uc_fast.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
uc_fast.set_conc({"A": 10., "B": 50., "C": 20.})

# %%

# %% [markdown]
# ### Run the reaction

# %%
uc_fast.single_compartment_react(initial_step=0.004, duration=0.06, variable_steps=True)

# %%
uc_fast.get_history()

# %% [markdown]
# ### Plots changes of concentrations with time

# %%
uc_fast.plot_history(show_intervals=True)

# %%

# %%

# %% [markdown]
# # 2. Let's now repeat the simulation with the "slow" preset, which yields more data points
# ### and, generally, more accuracy

# %%
# Instantiate the simulator and specify the chemicals
uc_slow = UniformCompartment(reactions=uc_fast.get_reactions(), preset="slow")
# Re-use the chemicals and reactions of part 1, but now with the "slow" preset

# %%
# Set the initial concentrations of all the chemicals
uc_slow.set_conc({"A": 10., "B": 50., "C": 20.})

# %%
uc_slow.single_compartment_react(initial_step=0.004, duration=0.06, variable_steps=True)

# %% [markdown]
# ### Note that **80 steps** were now used, instead of the earlier 23

# %%
uc_slow.get_history()

# %%

# %%

# %% [markdown]
# # 3. And, finally, get the EXACT analytical solution

# %%
# Instantiate the simulator and specify the chemicals
uc_exact = UniformCompartment(preset="slow", reactions=uc_fast.get_reactions(),
                              exact=True)
# Re-use the chemicals and reactions of part 1; 
# the "slow" preset doesn't particularly matter now because it's just 1 reaction, and we'll be getting the analytical solutin.  
# Note the `exact` flag

# %%
# Set the initial concentrations of all the chemicals
uc_exact.set_conc({"A": 10., "B": 50., "C": 20.}, snapshot=True)

# %%
uc_exact.single_compartment_react(initial_step=0.004, duration=0.06, variable_steps=True)

# %%
uc_exact.get_history() 

# %%

# %%

# %% [markdown]
# # 4. Let's compare plots for the concentration of `C`, as a function of time, from the earlier 3 simulations

# %%
p1 = uc_fast.plot_history(chemicals="C", colors=['#F5B914'], title="fast (less precise)", show=False)
p2 = uc_slow.plot_history(chemicals="C", colors=['#CBE504'], title="slow (more precise)", show=False)
p3 = uc_exact.plot_history(chemicals="C", colors=['forestgreen'], title="exact",  show=False)

# %%
PlotlyHelper.combine_plots([p1, p2, p3], title="Comparison of simulation accuracies", xrange=[0, 0.015])

# %% [markdown]
# ## The comparison reveals a gradient "less precise" -> "more precise" -> "exact" solution

# %%
