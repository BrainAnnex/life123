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
# ### IRREVERSIBLE unimolecular reaction `A -> B`,
# with 1st-order kinetics.
#
# **Adaptive variable time steps** 
# compared with **exact analytical solution**

# %% [markdown]
# ### TAGS :  "uniform compartment", "numerical", "quick-start", "basic", "under-the-hood"

# %%
LAST_REVISED = "Mar. 18, 2026"
LIFE123_VERSION = "1.0.0rc7"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import numpy as np

from life123 import check_version, UniformCompartment, ReactionKinetics, GraphicLog, PlotlyHelper

# %%
check_version(LIFE123_VERSION)    # To check compatibility

# %%

# %%

# %% [markdown]
# # PART 1 - VARIABLE TIME STEPS (Numerical Approximation to Solution)

# %%
# Instantiate the simulator and specify the chemicals
# Here we use the "fast" preset for the variable steps, 
# trying to push the envelope on speed
uc_approx = UniformCompartment(preset="fast")

# Reaction A <-> B , with 1st-order kinetics in both directions
uc_approx.add_reaction(reactants="A", products="B", 
                       kF=3., reversible=False) 

uc_approx.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
uc_approx.set_conc({"A": 50., "B": 10.})
uc_approx.describe_state()

# %%
uc_approx.get_history()

# %%
uc_approx.enable_diagnostics()   # To save diagnostic information about the call to single_compartment_react()
                          # Useful for insight into the inner workings of the simulation;
                          # in particular, for using plot_step_sizes()

# %%

# %% [markdown]
# ## Run the reaction   
# #### Passing True (default) to _variable_steps_ automatically adjusts up or down the time steps

# %%
uc_approx.single_compartment_react(initial_step=0.1, target_end_time=1.5,
                                   variable_steps=True)

# %%
history = uc_approx.get_history()   # The system's history, saved during the run of single_compartment_react()
history

# %%
uc_approx.plot_history(show_intervals=True)    
# We'll let the system pick default colors and assign them to the various chemicals

# %% [markdown]
# ### Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# That resulted from passing the flag _variable_steps=True_ to single_compartment_react()

# %%
uc_approx.plot_step_sizes(show_intervals=True)        #  To see the sizes of the steps taken

# %% [markdown]
# Why the zigzag?  It's because of the **"fast" preset** picked for the variable steps, in the instantiation of the class `UniformCompartment`: it's like a "high-strung driver" that tries to get away with excessive speeed - and periodically overdoes on acceleration, and then slams on the brakes!  
# Other presets (such as "mid") are more "mild-mannered" and more conservative about going too fast too soon.

# %%

# %%

# %% [markdown]
# # PART 2 - Exact analytical solution

# %%
uc_exact = UniformCompartment(reactions=uc_approx.get_reactions(), 
                              preset="fast", exact=True)   
# Re-use the chemicals and reactions of part 1 , and again employ the "fast" preset.  
# Note the `exact` flag

# %%
# Set the initial concentrations of all the chemicals
uc_exact.set_conc({"A": 50., "B": 10.})
uc_exact.describe_state()

# %%
uc_exact.enable_diagnostics() # To save diagnostic information about the call to single_compartment_react()
                              # Useful for insight into the inner workings of the simulation;
                              # in particular, for using plot_step_sizes()

# %% [markdown]
# ### Run the reaction (still with the default adaptive variable time steps, but now using the _exact_ solver)

# %%
uc_exact.single_compartment_react(initial_step=0.1, target_end_time=1.5,
                                  variable_steps=True)

# %%
uc_exact.plot_history(show_intervals=True, title_prefix="EXACT solution")

# %%
uc_exact.plot_step_sizes(show_intervals=True)        #  To see the sizes of the steps taken

# %% [markdown]
# A tad less "high-strung" than before...   
# Even though - with just 1 reaction and nothing else - all the values are absolutely exact, the simulator doesn't make special allowance for such a fact (of no use other than for a demo!) and still reins in situations that trigger various norms (measures on the rates of change, etc.)   
# Every simulation preset makes use of a variety of norms, for assessing adjustments to adaptive timestep sizes.

# %%

# %%

# %% [markdown]
# # PART 3 - Comparison of the two solutions

# %% [markdown]
# #### To avoid clutter, we'll just plot [A], as obtained from the approx solution and the exact analytical solution

# %%
fig_approx = uc_approx.plot_history(chemicals='A', colors='red', title="Approx solution", show=True)     
# Repeat a portion of the diagram seen in Part 1

# %%
fig_exact = uc_exact.plot_history(chemicals='A', colors='green', title="EXACT solution", show=True)     
# Repeat a portion of the diagram seen in Part 2

# %%
# Place both plots together

PlotlyHelper.combine_plots(fig_list=[fig_approx, fig_exact],
                           y_label="concentration [A]",
                           title="Approx vs. Exact soln, for [A] in irreversible reaction `A->B`",
                           legend_title="Simulation run")    

# %% [markdown]
# ### A pretty good overlap!

# %%
