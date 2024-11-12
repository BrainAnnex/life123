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
# ## Comparison of:  
# (1) adaptive variable time steps  
# (2) fixed time steps   
# (3) exact solution  
# ### for the reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium.
#
# This is a continuation of the experiments `react_2_a` (fixed time steps) and  `react_2_b` (adaptive variable time steps)
#
# **Background**: please see experiments `react_2_a` and `react_2_b`   

# %%
LAST_REVISED = "Nov. 11, 2024"
LIFE123_VERSION = "1.0.0.rc.0"      # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import check_version, ChemData, UniformCompartment

import numpy as np
import plotly.graph_objects as go
from life123 import ReactionKinetics
from life123.visualization.plotly_helper import PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%

# %%

# %% [markdown]
# ## Common set up for the chemicals and the reaction (used by all the simulations)

# %% tags=[]
# Instantiate the simulator and specify the chemicals
chem_data = ChemData()

# Reaction A <-> B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants="A", products="B", 
                       forward_rate=3., reverse_rate=2.)

chem_data.describe_reactions()

# %%

# %%

# %% [markdown]
# # PART 1 - VARIABLE TIME STEPS
# We'll do this part first, because the number of steps taken is unpredictable;  
# we'll note that number, and in Part 2 we'll do exactly that same number of fixed steps

# %%
dynamics_variable = UniformCompartment(chem_data=chem_data, preset="mid")

# Initial concentrations of all the chemicals
dynamics_variable.set_conc({"A": 10., "B": 50.})

dynamics_variable.describe_state()

# %% [markdown] tags=[]
# ### Run the reaction (VARIABLE adaptive time steps)

# %%
dynamics_variable.single_compartment_react(initial_step=0.1, target_end_time=1.2,
                                  variable_steps=True)

# %% [markdown]
# #### The flag _variable_steps_ automatically adjusts up or down the time steps
# In part 2, we'll remember that it took 19 steps

# %%
dynamics_variable.get_history()   # The system's history, saved during the run of single_compartment_react()

# %%
dynamics_variable.plot_history(colors=['darkturquoise', 'green'], show_intervals=True)

# %% [markdown] tags=[]
# #### Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# That resulted from passing the flag _variable_steps=True_ to single_compartment_react()

# %%

# %%

# %% [markdown]
# # PART 2 - FIXED TIME STEPS

# %% [markdown]
# #### This is a re-do of the above simulation simulation, but with a fixed time step
# The fixed time step is chosen to attain the same total number of data points as obtained with the variable time steps of part 1

# %%
dynamics_fixed = UniformCompartment(chem_data=chem_data)   # Re-use same chemicals and reactions of part 1

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
# # PART 3 - EXACT Solution

# %%
t_arr = np.linspace(0., 1.2, 41)   # A relatively dense uniform grid across our time range
t_arr

# %%
# The exact solution is available for a simple scenario like the one we're simulating here
rxn = chem_data.get_reaction(0)   # Object of type life123.reaction.Reaction

A_exact, B_exact = ReactionKinetics.solve_exactly(rxn=rxn, A0=10., B0=50., t_arr=t_arr)

# %%
fig_exact = PlotlyHelper.plot_curves(x=t_arr, y=[A_exact, B_exact], title="EXACT solution", x_label="SYSTEM TIME", y_label="concentration",
                                     legend_title="Chemical", curve_labels=["A (EXACT)", "B (EXACT)"],
                                     colors=["darkturquoise", "green"], show=True)

# %%

# %%

# %% [markdown] tags=[]
# # PART 4 - Comparing Variable Steps, Fixed Steps and Exact Solution   
# #### To avoid clutter, we'll just plot [A]

# %%
fig_variable = dynamics_variable.plot_history(chemicals=['A'], colors='darkturquoise', title="VARIABLE time steps", show=True)     # Repeat a portion of the diagram seen in Part 1

# %%
fig_fixed = dynamics_fixed.plot_history(chemicals=['A'], colors='blue', title="FIXED time steps", show=True)         # Repeat a portion of the diagram seen in Part 2

# %%
fig_exact = PlotlyHelper.plot_curves(x=t_arr, y=A_exact, title="EXACT solution", x_label="SYSTEM TIME", y_label="concentration",
                                     curve_labels="A (EXACT)", legend_title="Chemical",
                                     colors="red", show=True)  # Repeat a portion of the diagram seen in Part 3

# %%
PlotlyHelper.combine_plots(fig_list=[fig_fixed, fig_variable, fig_exact],
                           xrange=[0, 0.4], y_label="concentration [A]",
                           title="Variable time steps vs. Fixed vs. Exact soln, for [A] in reaction `A<->B`",
                           legend_title="Simulation run")    # All the 3 plots put together: show only the initial part (but it's all there; you can zoom out!)

# %% [markdown]
# #### Not surprisingly, the adaptive variable time steps outperform the fixed ones (for the same total number of points in the time grid), at times when there's pronounced change.  
# If you zoom out on the plot (by hovering on it, and using the Plotly controls that appear on the right, above), you can see all 3 curves essentially converging as the reaction approaches equilibrium.

# %%

# %%

# %% [markdown] tags=[]
# # PART 5 - Repeating Part 4 with a coarser grid
# #### The advantage of adaptive variable step will be even more prominent

# %%
# A coarser version of the variable-step simulation of Part 1
dynamics_variable_new = UniformCompartment(chem_data=chem_data, preset="fast")   # Re-use same chemicals and reactions of part 2

dynamics_variable_new.set_conc([10., 50.])

# %%
dynamics_variable_new.single_compartment_react(initial_step=0.1, target_end_time=1.2,
                                              variable_steps=True)

# %% [markdown]
# ### Note that the variable-step simulation is now taking 14 steps instead of the earlier 19

# %%
fig_variable = dynamics_variable_new.plot_history(chemicals='A', colors='darkturquoise', title="VARIABLE time steps",
                                              show_intervals=True, show=True)

# %%
# Now, a coarser version of the fixed-step simulation of Part 2
dynamics_fixed_new = UniformCompartment(chem_data=dynamics_fixed.chem_data)   # Re-using same chemicals and reactions of part 2

dynamics_fixed_new.set_conc([10., 50.])

# %%
# Matching the NEW total number of steps
dynamics_fixed_new.single_compartment_react(n_steps=14, target_end_time=1.2,
                                        variable_steps=False,
                                        snapshots={"initial_caption": "1st reaction step",
                                                   "final_caption": "last reaction step"})

# %%
fig_fixed = dynamics_fixed_new.plot_history(chemicals='A', colors='blue', title="FIXED time steps",
                                        show_intervals=True, show=True)

# %% [markdown]
# #### Notice the jaggedness at the left (jaggedness NOT present with the same number of total grid points, with the variable-step simulation)

# %%
PlotlyHelper.combine_plots(fig_list=[fig_fixed, fig_variable, fig_exact],
                           curve_labels = ["FIXED time steps", "VARIABLE time steps", "EXACT solution"],
                           xrange=[0, 0.4], y_label="concentration [A]",
                           title="Fixed vs. Variable time steps vs. Exact soln, for [A] in reaction `A<->B`",
                           legend_title="Simulation run")

# %% [markdown]
# ### With fewer grid points, the advantage of adaptive variable timesteps is more pronounced  
# If you zoom out the plot, and scroll to later times, you can see that the advantage later disappears when there's "less happening" (change-wise), closer to equilibrium

# %%
