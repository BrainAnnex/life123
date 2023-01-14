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
# ## A <-> B , shunted by coupled reaction A <-> S
# **Scenario 1** : The shunt (A <-> S) has a *kinetic* advantage but *thermodynamic* DISadvantage compared to A <-> B   
# (i.e. A <-> S is fast, but energetically unfavored) 
#
# **Scenario 2** : The shunt (A <-> S) is has a *kinetic* DISadvantage but a *thermodynamic* advantage compared to A <-> B     
# (i.e. A <-> S is slow, but energetically favored)  
#
# All reactions 1st order, mostly forward.  Taken to equilibrium.
#
# LAST REVISED: Jan. 13, 2023

# %% [markdown]
# # TODO: maybe also show A <-> B in the absence of the 2nd reaction

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

# Reaction A <-> B
chem_data.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=30., reverse_rate=5.) 

# Reaction A <-> S (fast shunt, poor thermodynical energetic advantage)
chem_data.add_reaction(reactants=["A"], products=["S"],
                       forward_rate=150., reverse_rate=100.) 

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
#dynamics.verbose_list = [1, 2, 3]      # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; 
# so, we'll be using the dynamic_steps option to increase time resolution,
# as long as the reaction remains "fast" (based on a threshold of % change, as specified by fast_threshold)
dynamics.single_compartment_react(time_step=0.001, reaction_duration=0.05,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_steps=4, fast_threshold=10)

# %%
# Continue running the reaction at lover resolution
dynamics.single_compartment_react(time_step=0.002, reaction_duration=0.25,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_steps=4, fast_threshold=10)

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
df = dynamics.history.get()
df

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=3)

# %% [markdown] tags=[]
# ## <a name="react_5_plot"></a>Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "B", "S"], 
              title="Coupled reactions A <-> B and A <-> S",
              color_discrete_sequence = ['blue', 'green', 'red'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%

# %%
# Verify that the stoichiometry is respected at every reaction step/substep (NOTE: it requires earlier activation of saving diagnostic data)
dynamics.stoichiometry_checker_entire_run()

# %%

# %%

# %% [markdown]
# # Re-Start the simulation
# ## slow shunt, with excellent thermodynamical energetic advantage

# %%
# Specify the chemicals
chem_data2 = chem(names=["A", "B", "S"])

# Reaction A <-> B
chem_data2.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=30., reverse_rate=5.) 

# Reaction A <-> S (slow shunt, excellent thermodynamical energetic advantage)
chem_data2.add_reaction(reactants=["A"], products=["S"],
                       forward_rate=3., reverse_rate=0.1)

chem_data2.describe_reactions()

# %%
dynamics2 = ReactionDynamics(reaction_data=chem_data2)

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics2.set_conc([50., 0, 0.], snapshot=True)

# %%
dynamics2.describe_state()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics2.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()
#dynamics.verbose_list = [1, 2, 3]      # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; 
# so, we'll be using the dynamic_steps option to increase time resolution,
# as long as the reaction remains "fast" (based on a threshold of % change, as specified by fast_threshold)
dynamics2.single_compartment_react(time_step=0.005, reaction_duration=0.3,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_steps=5, fast_threshold=10)

# %%
# Continue running the reaction at lover resolution
dynamics2.single_compartment_react(time_step=0.25, reaction_duration=6.7,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_steps=5, fast_threshold=10)

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
df2 = dynamics2.history.get()
df2

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics2.is_in_equilibrium(tolerance=13)

# %% [markdown] tags=[]
# ## <a name="react_5_plot"></a>Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics2.get_history(), x="SYSTEM TIME", y=["A", "B", "S"], 
              title="Coupled reactions A <-> B and A <-> S",
              color_discrete_sequence = ['blue', 'green', 'red'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
