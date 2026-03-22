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
# ## Comparison of:  
# (1) adaptive variable time steps  
# (2) fixed time steps   
# (3) exact solution  
# ### for the elementary reaction `A <-> B`,
# taken to equilibrium.
#
# This is a continuation of the experiments `react_2_a` (fixed time steps) and  `react_2_b` (adaptive variable time steps)
#
# **Background**: please see experiments `react_2_a` and `react_2_b`   

# %% [markdown]
# ### TAGS : "numerical", "uniform compartment"

# %%
LAST_REVISED = "Mar. 18, 2026"
LIFE123_VERSION = "1.0.0rc7"        # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import check_version, ChemData, UniformCompartment

import numpy as np
from life123 import ReactionKinetics, ReactionRegistry, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%

# %%

# %% [markdown]
# ## Common set up for the chemicals and the reaction (used by all the simulations)

# %%
# Set up the reactions and their chemicals (common for all the simulations below)
rxns = ReactionRegistry()

rxns.add_reaction(reactants="A", products="B", kF=3., kR=2.)  # ELementary reaction A <-> B

rxns.describe_reactions()

# %%

# %% [markdown]
# # PART 1 - VARIABLE TIME STEPS
# We'll do this part first, because the number of steps taken is unpredictable;  
# we'll note that number, and in Part 2 we'll do exactly that same number of fixed steps

# %%
dynamics_variable = UniformCompartment(reactions=rxns, preset="mid")

# Initial concentrations of all the chemicals
dynamics_variable.set_conc({"A": 10., "B": 50.})

dynamics_variable.describe_state()

# %% [markdown]
# ### Run the reaction (VARIABLE adaptive time steps)

# %%
dynamics_variable.single_compartment_react(initial_step=0.1, target_end_time=1.2,
                                           variable_steps=True)
# Specifying variable_steps=True for emphasis; it's the default

# %% [markdown]
# #### The flag _variable_steps_ automatically adjusts up or down the time steps
# In part 2, we'll remember that it took **19 steps**

# %%
dynamics_variable.get_history()   # The system's history, saved during the run of single_compartment_react()

# %%
dynamics_variable.plot_history(colors=['darkturquoise', 'green'], show_intervals=True)

# %% [markdown]
# #### Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# That resulted from passing the flag _variable_steps=True_ to  `single_compartment_react()`

# %%

# %%

# %% [markdown]
# # PART 2 - FIXED TIME STEPS

# %% [markdown]
# #### This is a re-do of the above simulation, but with a fixed time step
# The fixed time step is chosen to attain the same total number of data points (i.e. 19) as obtained with the variable time steps of part 1

# %%
dynamics_fixed = UniformCompartment(reactions=rxns)   # Re-use the chemicals and reactions of part 1

# %%
# Initial concentrations of all the chemicals
dynamics_fixed.set_conc({"A": 10., "B": 50.})

dynamics_fixed.describe_state()

# %% [markdown]
# ### Run the reaction (FIXED time steps)

# %%
# Matching the total number of steps to the earlier, variable-step simulation
dynamics_fixed.single_compartment_react(n_steps=19, target_end_time=1.2,
                                        variable_steps=False)

# %%
dynamics_fixed.get_history()   # The system's history, saved during the run of single_compartment_react()

# %%
dynamics_fixed.plot_history(colors=['darkturquoise', 'green'], show_intervals=True)

# %% [markdown]
# Notice how grid points are being "wasted" on the tail part of the simulation, where little is happening - grid points that would be best used in the early part, as was done by the variable-step simulation of Part 1

# %%

# %%

# %% [markdown]
# # PART 3 - EXACT Solution (same fixed steps as in Part 2)

# %%
dynamics_exact = UniformCompartment(reactions=rxns, exact=True)   
# Re-use the chemicals and reactions of part 1 .  Note the `exact` flag

# %%
# Initial concentrations of all the chemicals
dynamics_exact.set_conc({"A": 10., "B": 50.})

dynamics_exact.describe_state()

# %% [markdown]
# ### Run the reaction (same FIXED time steps as before, but now using the exact solver)

# %%
# Matching the total number of steps to those of Part 1 and Part 2
dynamics_exact.single_compartment_react(n_steps=19, target_end_time=1.2,
                                        variable_steps=False)

# %%
dynamics_exact.get_history()   # The system's history, saved during the run of single_compartment_react()

# %%
dynamics_exact.plot_history(colors=['darkturquoise', 'green'], show_intervals=True, title_prefix="EXACT solution")

# %%

# %%

# %% [markdown]
# # PART 4 - Comparing Variable Steps, Fixed Steps and Exact Solution   
# #### To avoid clutter, we'll just plot [A]

# %%
# A streamlined version of the diagram seen in Part 1
fig_variable = dynamics_variable.plot_history(chemicals='A', colors='darkturquoise', title="VARIABLE time steps", show=True)     

# %%
# A streamlined version of the diagram seen in Part 2
fig_fixed = dynamics_fixed.plot_history(chemicals='A', colors='blue', title="FIXED time steps", show=True)

# %%
# A streamlined version of the diagram seen in Part 3
fig_exact = dynamics_exact.plot_history(chemicals='A', colors='red', title="EXACT solution (at fixed time steps)", show=True)

# %%
PlotlyHelper.combine_plots(fig_list=[fig_fixed, fig_variable, fig_exact],
                           xrange=[0, 0.8], y_label="concentration [A]",
                           title="Variable time steps vs. Fixed vs. Exact soln, for [A] in reaction `A<->B`",
                           legend_title="Simulation run")    
# All the 3 plots put together: show only the initial part (but it's all there; you can zoom out!)

# %% [markdown]
# #### Not surprisingly, the adaptive variable time steps outperform the fixed ones (for the same total number of points in the time grid), in the early/middle part of the plot - i.e. at times when there's pronounced change.  
# All 3 curves essentially converge as the reaction approaches equilibrium.

# %%
