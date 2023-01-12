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
# ## 2 COUPLED reactions of different speeds:  
# ### A <-> B (fast) and A <-> S (slow)
# All 1st order. Taken to equilibrium. Both reactions are mostly forward.
#
# (Adaptive variable time resolution is used, with extensive diagnostics.)
#
# LAST REVISED: Jan. 10, 2023

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from modules.reactions.reaction_data import ReactionData as chem
from modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
from modules.visualization.graphic_log import GraphicLog

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
chem_data = chem(names=["A", "B", "S"])

# Reaction A <-> B (fast)
chem_data.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=100., reverse_rate=8.) 

# Reaction A <-> S (slow)
chem_data.add_reaction(reactants=["A"], products=["S"],
                       forward_rate=6., reverse_rate=2.) 

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

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics.set_conc([50., 0, 0.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()
dynamics.verbose_list = [1, 2, 3]      # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; 
# so, we'll be using the dynamic_steps option to increase time resolution,
# as long as the reaction remains "fast" (based on a threshold of % change, as specified by fast_threshold)
dynamics.single_compartment_react(time_step=0.02, reaction_duration=0.8,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_steps=10, fast_threshold=20)

# %%
df = dynamics.history.get()
df

# %%
df.loc[:50]

# %%
df.loc[51:99]

# %%
# Let's expand the last part
df.loc[100:]

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=0.4)

# %% [markdown] tags=[]
# ## <a name="react_5_plot"></a>Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "B", "S"], 
              title="Coupled reactions A <-> B and A <-> S",
              color_discrete_sequence = ['blue', 'red', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%

# %%

# %% [markdown]
# # EVERYTHING BELOW IS FOR DIAGNOSTIC INSIGHT
#
# ### Perform some verification

# %%
# Verify that the stoichiometry is respected at every reaction step/substep (NOTE: it requires earlier activation of saving diagnostic data)
dynamics.stoichiometry_checker_entire_run()

# %% [markdown]
# ### Take a peek at the diagnostic data saved during the earlier reaction simulation

# %%
# This table contains a subset of a typical system history
dynamics.diagnostic_data_baselines.get()

# %%
# Concentration increments due to reaction 0 (A <-> B)
# Note that [S] is not affected
dynamics.get_diagnostic_data(rxn_index=0)

# %%
# Expand the last part of the above table
dynamics.get_diagnostic_data(rxn_index=0).loc[60:]

# %%
# Concentration increments due to reaction 1 (A <-> S)
# Note that [B] is not affected
dynamics.get_diagnostic_data(rxn_index=1)

# %%

# %% [markdown]
# ### Provide a detailed explanation of all the steps/substeps of the reactions, from the saved diagnostic data

# %%
dynamics.explain_reactions()

# %%

# %%

# %%

# %% [markdown]
# # Re-run with constanst steps

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)

# %%
dynamics.set_conc([50., 0, 0.], snapshot=True)

# %%
dynamics.single_compartment_react(time_step=0.005, reaction_duration=0.8,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  )      

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "B", "S"], 
              title="Coupled reactions A <-> B and A <-> S",
              color_discrete_sequence = ['blue', 'red', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=0.4)

# %%

# %%
df = dynamics.history.get()
df

# %%
