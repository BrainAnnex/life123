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
# LAST REVISED: Nov. 26, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import plotly.graph_objects as go

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
# ### This is a re-do of the simulation of `react_1_a` but with a smaller time step
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
history_fixed = dynamics_fixed.get_history()   # The system's history, saved during the run of single_compartment_react()
history_fixed

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
history_variable = dynamics_variable.get_history()   # The system's history, saved during the run of single_compartment_react()
history_variable

# %%
dynamics_variable.plot_history(colors=['blue', 'orange'], show_intervals=True)

# %% [markdown] tags=[]
# #### Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# That resulted from passing the flag _variable_steps=True_ to single_compartment_react()

# %%

# %% [markdown] tags=[]
# # PART 3 - Comparing Fixed Steps, Variable Steps and Exact Solution

# %%
A_fixed = history_fixed["A"].to_numpy()
A_fixed

# %%
A_variable = history_variable["A"].to_numpy()
A_variable

# %%
len(A_fixed)

# %%
len(A_variable)

# %%
fig_fixed = dynamics_fixed.plot_history(colors=['blue', 'orange'], title="FIXED time steps")
fig_fixed

# %%
fig_variable = dynamics_variable.plot_history(colors=['navy', 'brown'], title="VARIABLE time steps")
fig_variable

# %%
all_fig = go.Figure(data=fig_fixed.data + fig_variable.data, 
                    layout = fig_fixed.layout)    # Note that the + is concatenating lists
all_fig.update_layout(title="Fixed vs. Variable time steps, for reaction `A<->B`")

# Set the names to show in the legend
all_fig['data'][0]['name']='A (fixed)'
all_fig['data'][1]['name']='B (fixed)' 
all_fig['data'][2]['name']='A (variable)' 
all_fig['data'][3]['name']='B (variable)' 

all_fig.show()

# %% [raw]
#

# %%
