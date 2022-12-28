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
# ## Association/Dissociation reaction 2A <-> C
# #### with 2nd-order kinetics for A,  
# #### and 1-st order kinetics for C
#
# Compare with experiment "react_3"
#
#
# LAST REVISED: Dec. 25, 2022  ************* IN-PROGRESS  **********

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
chem_data = chem(names=["A", "C"])

# Reaction 2A <-> C , with 2nd-order kinetics for A, and 1st-order kinetics for C
chem_data.add_reaction(reactants=[(2, "A", 2)], products=["C"],
                       forward_rate=5., reverse_rate=2.)   
# Note: the first 2 in (2, "A", 2) is the stoichiometry coefficient, while the other one is the order

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
dynamics.set_conc([60., 20.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.single_compartment_react(time_step=0.001, reaction_duration=0.06,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_step=2)      
                                  # Accepting the default:  fast_threshold=5

# %% [markdown]
# ### Note: the argument _dynamic_step=2_ splits the time steps in 2 whenever the reaction is "fast" (as determined using fast_threshold=5)

# %%
df = dynamics.history.get()
df

# %% [markdown]
# ### Notice how the reaction proceeds in smaller steps in the early times, when the concentrations are changing much more rapidly

# %% [markdown]
# ## Note: "A" (now largely depleted) is largely the limiting reagent

# %% [markdown]
# ### Check the final equilibrium

# %%
dynamics.get_system_conc()

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(rxn_index=0, conc=dynamics.get_conc_dict())

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "C"], 
              title="Reaction 2A <-> C  (2nd order in A).  Changes in concentrations with time",
              color_discrete_sequence = ['red', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
