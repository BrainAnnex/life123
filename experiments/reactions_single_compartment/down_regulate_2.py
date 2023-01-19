# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# ## A down-regulates B , 
# ### by being the *limiting reagent* in reaction A + 2 B <-> Y (mostly forward)
# 1st-order kinetics.   
# If [A] is low, [B] remains high.  Then, if [A] goes high, [B] goes low.  However, at that point, A can no longer bring B up to any substantial extent.
#
# See also 1D/reactions/down_regulation_1
#
# LAST REVISED: Jan. 18, 2023

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
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# ### Initialize the system

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "Y"])

# Reaction A + X <-> 2B , with 1st-order kinetics for all species
chem_data.add_reaction(reactants=[("A") , (2, "B")], products=[("Y")],
                       forward_rate=8., reverse_rate=2.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc([5, 100., 0.], snapshot=True)  # A is scarce, B is plentiful, C is absent
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Take the initial system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.0005, n_steps=30)

df = dynamics.history.get()
df

# %% [markdown]
# A, as the scarse limiting reagent, stops the reaction.  
# When A is low, B is also low.

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['red', 'green', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown] tags=[]
# ## Now, let's suddenly increase [A]

# %%
dynamics.set_chem_conc(species_index=0, conc=40., snapshot=True)
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.0005, n_steps=80)

df = dynamics.history.get()
df

# %% [markdown]
# **A**, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=7)

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['red', 'green', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown] tags=[]
# ## Let's again suddenly increase [A]

# %%
dynamics.set_chem_conc(species_index=0, conc=30., snapshot=True)
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.0005, n_steps=70)

df = dynamics.history.get()
df

# %% [markdown]
# **A**, again the scarse limiting reagent, stops the reaction yet again

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['red', 'green', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# Note: A can up-regulate B, but it cannot bring it down.  
# X will soon need to be replenished, if A is to continue being the limiting reagent.

# %%
