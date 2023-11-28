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
# ### Comparison of fixed time steps, adaptive variable time steps, and exact solution 
# ### for reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# This is a continuation of the experiments _"react_1_a"_ (fixed time steps) and  _"react_1_b"_ (adaptive variable time steps)
#
# LAST REVISED: Nov. 27, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.graph_objects as go
from src.modules.visualization.plotly_helper import PlotlyHelper

# %%

# %% [markdown]
# ### Common set up of the chemicals and the reaction (used by all the simulations)

# %% tags=[]
# Specify the chemicals
chem_data = chem(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants=["A"], products=["B"], 
                       forward_rate=3., reverse_rate=2.)

chem_data.describe_reactions()

# %%

# %% [markdown]
# # PART 1 - FIXED TIME STEPS

# %% [markdown]
# #### This is a re-do of the simulation of `react_1_a` but with a smaller time step
# The fixed time step is chosen to attain the same total number of data points as the variable time steps of part 2

# %%
dynamics_fixed = ReactionDynamics(chem_data=chem_data)

# %%
# Initial concentrations of all the chemicals, in index order
dynamics_fixed.set_conc([10., 50.])

dynamics_fixed.describe_state()

# %% [markdown]
# ### Run the reaction (FIXED time steps)

# %%
dynamics_fixed.single_compartment_react(initial_step=0.06, target_end_time=1.2,
                                        variable_steps=False,
                                        snapshots={"initial_caption": "1st reaction step",
                                                   "final_caption": "last reaction step"})

# %%
dynamics_fixed.get_history()   # The system's history, saved during the run of single_compartment_react()

# %%
dynamics_fixed.plot_history(colors=['blue', 'orange'], show_intervals=True)

# %%

# %% [markdown]
# # PART 2 - VARIABLE TIME STEPS

# %%
dynamics_variable = ReactionDynamics(chem_data=chem_data)

# %%
# Initial concentrations of all the chemicals, in index order
dynamics_variable.set_conc([10., 50.])

dynamics_variable.describe_state()

# %%
# For experiment repeatability, we specify a particular group of preset parameters applicable to the adaptive time steps
dynamics_variable.set_adaptive_parameters(preset="mid")   # A "middle-of-the road" heuristic: somewhat "conservative" but not overly so

# %% [markdown] tags=[]
# ### Run the reaction (VARIABLE adaptive time steps)

# %%
dynamics_variable.single_compartment_react(initial_step=0.1, target_end_time=1.2,
                                  variable_steps=True, explain_variable_steps=False,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"}
                                  )

# %% [markdown]
# #### The flag _variable_steps_ automatically adjusts up or down the time step,  whenever the changes of concentrations are, respectively, "slow" or "fast" (as determined using the specified _thresholds_ )

# %%
dynamics_variable.get_history()   # The system's history, saved during the run of single_compartment_react()

# %%
dynamics_variable.plot_history(colors=['blue', 'orange'], show_intervals=True)

# %% [markdown] tags=[]
# #### Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# That resulted from passing the flag _variable_steps=True_ to single_compartment_react()

# %%

# %% [markdown] tags=[]
# # PART 3 - Comparing Fixed Steps, Variable Steps and Exact Solution

# %%
fig_fixed = dynamics_fixed.plot_history(chemicals=['A'], colors='blue', title="FIXED time steps", range_x=[0, 0.4])
fig_fixed

# %%
fig_variable = dynamics_variable.plot_history(chemicals=['A'], colors='aqua', title="VARIABLE time steps", range_x=[0, 0.4])
fig_variable

# %%
fig_2 = PlotlyHelper.combine_plots(fig_list=[fig_fixed, fig_variable],
                           xrange=[0, 0.4], ylabel="concentration [A]",
                           title="Fixed vs. Variable time steps, for reaction `A<->B`",
                           legend_title="Simulation run", show=True)

# %%
t_arr = np.linspace(0., 1.0, 41)   # A uniform grid across our time range
t_arr

# %%
A_exact, B_exact = dynamics_variable.solve_exactly(rxn_index=0, A0=10., B0=50., t_arr=t_arr)

# %%
fig_exact = PlotlyHelper.plot_curves(x=t_arr, y=A_exact, title="EXACT solution", xlabel="SYSTEM TIME", ylabel="concentration", 
                                     curve_labels="A (EXACT)", legend_title="Chemical",
                                     colors="red", show=True)

# %%
PlotlyHelper.combine_plots(fig_list=[fig_fixed, fig_variable, fig_exact],
                           xrange=[0, 0.4], ylabel="concentration [A]",
                           title="Fixed vs. Variable time steps vs. Exact soln, for [A] in reaction `A<->B`",
                           legend_title="Simulation run")

# %% [markdown]
# Not surprisingly, the adaptive variable time steps outperform the fixed ones (for an equal total number of grid points), at times when there's pronounced change.  
# If you zoom out on the plot (hover and use the controls on the right, above), you can see all 3 curves essentially converging as the reaction approaches equilibrium.

# %%
