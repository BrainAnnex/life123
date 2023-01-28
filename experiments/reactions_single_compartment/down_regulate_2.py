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
# ## 'A' down-regulates 'B' , 
# ### by being the *limiting reagent* in reaction A + 2 B <-> Y (mostly forward)
# 1st-order kinetics.   
# If [A] is low and [B] is high, then [B] remains high.  If [A] goes high, [B] goes low.  However, at that point, A can no longer bring B up to any substantial extent.
#
# See also 1D/reactions/down_regulation_1
#
# LAST REVISED: Jan. 28, 2023

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

# Reaction A + 2 B <-> Y , with 1st-order kinetics for all species
chem_data.add_reaction(reactants=[("A") , (2, "B")], products=[("Y")],
                       forward_rate=8., reverse_rate=2.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc={"A": 5., "B": 100., "Y": 0.},
                  snapshot=True)      # A is scarce, B is plentiful, Y is absent
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Take the initial system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.0005, n_steps=30, 
                                  dynamic_steps=2, fast_threshold=15)

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
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['red', 'blue', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown] tags=[]
# # Now, let's suddenly increase [A]

# %%
dynamics.set_chem_conc(species_name="A", conc=40., snapshot=True)
dynamics.describe_state()

# %%
dynamics.history.get(tail=5)

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.001, n_steps=40,
                                 dynamic_steps=2, fast_threshold=15)

df = dynamics.history.get()
df

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=7)

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations (reaction A + 2 B <-> Y)",
              color_discrete_sequence = ['red', 'blue', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# **A**, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %% [markdown] tags=[]
# # Let's again suddenly increase [A]

# %%
dynamics.set_chem_conc(species_name="A", conc=30., snapshot=True)
dynamics.describe_state()

# %%
dynamics.history.get(tail=5)

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.001, n_steps=35,
                                 dynamic_steps=2, fast_threshold=15)

df = dynamics.history.get()
df

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations (reaction A + 2 B <-> Y)",
              color_discrete_sequence = ['red', 'blue', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# **A**, again the scarse limiting reagent, stops the reaction yet again

# %%

# %% [markdown]
# ## A can down-regulate B, but it cannot bring it up to any significant amount
# #### Even if A is completely taken out (i.e., [A] set to 0), [B] can only slightly increase, from the reverse reaction ("Le Chatelier's principle".)  
# Let's try it:

# %%
dynamics.set_chem_conc(species_name="A", conc=0., snapshot=True)   # Completely eliminate A
dynamics.describe_state()

# %%
dynamics.single_compartment_react(time_step=0.001, n_steps=70,
                                 dynamic_steps=2, fast_threshold=15)

df = dynamics.history.get()
df

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations (reaction A + 2 B <-> Y)",
              color_discrete_sequence = ['red', 'blue', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# #### As expected, even the complete withdrawal of A (red), brings about only a modest increase of B's concentration, from the reverse reaction (i.e. [B] slightly increases at the expense of [Y].)  
# #### The change is modest because our  reaction A + 2 B <-> Y is mostly in the forward direction (K = 4)
# Le Chatelier's principle in action: "A change in one of the variables that describe a system at equilibrium produces a shift in the position of the equilibrium that counteracts the effect of this change."

# %%
