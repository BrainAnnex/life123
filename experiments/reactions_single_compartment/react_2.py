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
# ### Adaptive time steps (variable time resolution) for reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# Same as the experiment _"react_1"_ , but with adaptive variable time steps
#
# LAST REVISED: Mar. 12, 2023

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import pandas as pd
import plotly.express as px
from src.modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # Initialize the System
# Specify the chemicals and the reactions

# %% tags=[]
# Specify the chemicals
chem_data = chem(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants=["A"], products=["B"], 
                       forward_rate=3., reverse_rate=2.)

print("Number of reactions: ", chem_data.number_of_reactions())

# %%
chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# # Start the simulation

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)

# %%
# Initial concentrations of all the chemicals, in index order
dynamics.set_conc([10., 50.])

dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(initial_step=0.1, target_end_time=1.2,
                                  variable_steps=True, thresholds={"low": 0.5, "high": 0.8},
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  )
# Note: we are accepting the default variable_steps_threshold_abort = 1.44 and abort_step_reduction_factor = 2.

# %% [markdown]
# ## The flag _variable_steps_ automatically adjusts up or down the time step,  whenever the changes of concentrations are, respectively, "slow" or "fast" (as determined using the specified _thresholds_ )

# %%
df = dynamics.get_history()   # The system's history, saved during the run of single_compartment_react()
df

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# ### That resulted from passing the flag _variable_steps_ to single_compartment_react()

# %% [markdown] tags=[]
# ## Detailed Example 1: **going from 0.1375 to 0.1875**    

# %%
lookup = pd.merge_asof(pd.DataFrame({'lookup':[0.1375, 0.1875]}), df,
              right_on='SYSTEM TIME', left_on='lookup',
              direction='nearest')
lookup

# %%
delta_concentrations = (lookup.iloc[1][['A', 'B']] - 
                        lookup.iloc[0][['A', 'B']]).to_numpy()
delta_concentrations

# %% [markdown]
# As expected by the 1:1 stoichiometry, delta_A = - delta_B
#
# The above values coud also be looked up from the diagnostic data, since we only have 1 reaction:

# %%
rxn_data = dynamics.get_diagnostic_rxn_data(rxn_index=0)

# %%
rxn_data[0:12]

# %%
delta_row = dynamics.get_diagnostic_rxn_data(rxn_index=0, t=0.1375) # Locate the row with interval's start time
delta_row

# %%
delta_row[["Delta A", "Delta B"]].to_numpy()   # Gives same value as delta_concentrations, above

# %%
adjusted_L2_rate = dynamics.norm_A(delta_concentrations)  # A measure of how large delta_concentrations is
adjusted_L2_rate

# %%
dynamics.step_determiner_A(adjusted_L2_rate)

# %% [markdown]
# #### The above conclusion is that the step will be **HALVED** at the next round : that's because the adjusted_L2_rate > the "high" value given in the argument _thresholds={"low": 0.5, "high": 0.8}_

# %% [markdown] tags=[]
# ## Detailed Example 2: **going from 0.1875 to 0.2125**   

# %%
lookup = pd.merge_asof(pd.DataFrame({'lookup':[0.1875, 0.2125]}), df,
              right_on='SYSTEM TIME', left_on='lookup',
              direction='nearest')
lookup

# %%
delta_concentrations = (lookup.iloc[1][['A', 'B']] - 
                        lookup.iloc[0][['A', 'B']]).to_numpy()
delta_concentrations

# %% [markdown]
# Note how substantially smaller _delta_concentrations_ is, compared to the previous example

# %%
adjusted_L2_rate = dynamics.norm_A(delta_concentrations)  # A measure of how large delta_concentrations is
adjusted_L2_rate

# %%
dynamics.step_determiner_A(adjusted_L2_rate)

# %% [markdown]
# #### The above conclusion is that the step will be **DOUBLED** at the next round : that's because the adjusted_L2_rate < the "low" value given in the argument _thresholds={"low": 0.5, "high": 0.8}_

# %% [markdown]
# # Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['blue', 'orange'])

# %% [markdown]
# ## Note how the left-hand side of this plot is much smoother than it was in experiment `react_1`, where no adaptive time steps were used!

# %%
dynamics.plot_curves(colors=['blue', 'orange'], show_intervals=True)

# %%
dynamics.plot_step_sizes(show_intervals=True)

# %% [markdown]
# # Diagnostics of the run may be investigated as follows:  
# _(note - this is possible because we make a call to set_diagnostics() prior to running the simulation)_

# %%
dynamics.get_diagnostic_conc_data()

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)      # For the 0-th reaction (the only reaction in our case)

# %%
