# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
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
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "quick-start", "basic", "under-the-hood"

# %%
LAST_REVISED = "Oct. 6, 2024"
LIFE123_VERSION = "1.0.0.beta.39"   # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import ipynbname
import numpy as np

from life123 import check_version, UniformCompartment, ReactionKinetics, GraphicLog, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file
                                            # IN CASE OF PROBLEMS, set manually to any desired name

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%

# %%

# %% [markdown]
# # PART 1 - VARIABLE TIME STEPS (Numerical Approximation to Solution)

# %% tags=[]
# Instantiate the simulator and specify the chemicals
# Here we use the "fast" preset for the variable steps, trying to push the envelope on speed
uc = UniformCompartment(preset="fast")

# Reaction A <-> B , with 1st-order kinetics in both directions
uc.add_reaction(reactants="A", products="B", 
                forward_rate=3., reverse_rate=0)      # Notice the zero reverse rate (irreversible)

uc.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")

# %%
# Set the initial concentrations of all the chemicals
uc.set_conc({"A": 50., "B": 10.})
uc.describe_state()

# %%
uc.get_history()

# %%
uc.enable_diagnostics()   # To save diagnostic information about the call to single_compartment_react()
                          # Useful for insight into the inner workings of the simulation

# %% [markdown] tags=[]
# ## Run the reaction   
# #### Passing True (default) to _variable_steps_ automatically adjusts up or down the time steps

# %%
uc.single_compartment_react(initial_step=0.1, target_end_time=1.5,
                            variable_steps=True)

# %%
history = uc.get_history()   # The system's history, saved during the run of single_compartment_react()
history

# %%
uc.plot_history(colors=['darkturquoise', 'orange'], show_intervals=True)    # Plots of concentration with time

# %% [markdown] tags=[]
# ### Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# That resulted from passing the flag _variable_steps=True_ to single_compartment_react()

# %%
uc.plot_step_sizes(show_intervals=True)        #  To see the sizes of the steps taken

# %% [markdown]
# Why the zigzag?  It's because of the **"fast" preset** picked for the variable steps, in the instantiation of the class `UniformCompartment`: it's like a "high-strung driver" that tries to get away with excessive speeed - and periodically overdoes on acceleration, and then slams on the brakes!  
# Other presets (such as "mid") are more "mild-mannered" and more conservative about going too fast too soon.

# %%

# %%

# %% [markdown]
# # PART 2 - Comparison with exact analytical solution

# %%
t_arr = np.linspace(0., 1.5, 50)   # A relatively dense uniform grid across our time range
t_arr

# %%
# The exact solution is available for a simple scenario like the one we're simulating here

A_exact, B_exact = ReactionKinetics.exact_solution_unimolecular_irreversible(kF=3., A0=50., B0=10., t_arr=t_arr)

# %%
fig_exact = PlotlyHelper.plot_curves(x=t_arr, y=[A_exact, B_exact], title="EXACT solution", x_label="SYSTEM TIME", y_label="concentration",
                                     legend_title="Chemical", curve_labels=["A (EXACT)", "B (EXACT)"],
                                     colors=["darkturquoise", "orange"], show=True)

# %% [markdown] tags=[]
# #### To avoid clutter, we'll just plot [A], as obtained from the variable-step approx solution and the exact analytical solution

# %%
fig_exact = PlotlyHelper.plot_curves(x=t_arr, y=A_exact, title="EXACT solution", x_label="SYSTEM TIME", y_label="concentration",
                                     curve_labels="A (EXACT)", legend_title="Chemical",
                                     colors="red", show=True)     # Repeat a portion of the diagram seen just before

# %%
fig_variable = uc.plot_history(chemicals=['A'], colors='darkturquoise', title="VARIABLE time steps", show=True)     # Repeat a portion of the diagram seen in Part 1

# %%
PlotlyHelper.combine_plots(fig_list=[fig_variable, fig_exact],
                           xrange=[0, 1.5], ylabel="concentration [A]",
                           title="Variable time steps vs. Exact soln, for [A] in irreversible reaction `A->B`",
                           legend_title="Simulation run")    # Both plots put together

# %% [markdown]
# ### A pretty good overlap!

# %%
