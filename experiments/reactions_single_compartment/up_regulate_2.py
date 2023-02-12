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
# ## `U` ("Up-regulator") up-regulates `X` , by sharing a reaction product `D` ("Drain") across 2 separate reactions:   
# ### `U <-> 2 D` and `X <-> D` (both mostly forward)
#
# Invoking [Le Chatelier's principle](https://www.chemguide.co.uk/physical/equilibria/lechatelier.html), it can be seen that, starting from equilibrium, when [U] goes up, so does [D]; and when [D] goes up, so does [X].   
# Conversely, when [U] goes down, so does [D]; and when [D] goes down, so does [X].   
#
# 1st-order kinetics throughout.   
#
# LAST REVISED: Feb. 11, 2023

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
import plotly.express as px
from src.modules.visualization.graphic_log import GraphicLog

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
chem_data = chem(names=["U", "X", "D"])

# Reaction U <-> 2D , with 1st-order kinetics for all species
chem_data.add_reaction(reactants="U", products=[(2, "D")],
                       forward_rate=8., reverse_rate=2.)

# Reaction X <-> D , with 1st-order kinetics for all species
chem_data.add_reaction(reactants="X", products="D",
                       forward_rate=6., reverse_rate=3.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "D": 0.})
dynamics.describe_state()

# %% [markdown] tags=[]
# # 1. Take the initial system to equilibrium

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.03, stop_time=0.5,
                                  dynamic_substeps=2, rel_fast_threshold=50)

df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['red', 'green', 'gray'])

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# # 2. Now, let's suddenly increase [U]

# %%
dynamics.set_chem_conc(species_name="U", conc=70., snapshot=True)
dynamics.describe_state()

# %%
dynamics.get_history(tail=3)

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.03, stop_time=1,
                                  dynamic_substeps=2, rel_fast_threshold=50)

df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['red', 'green', 'gray'])

# %% [markdown]
# ### The (transiently) high value of [U] led to an increase in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

# %% [markdown] tags=[]
# # 3. Let's again suddenly increase [U]

# %%
dynamics.set_chem_conc(species_name="U", conc=100., snapshot=True)
dynamics.describe_state()

# %%
dynamics.get_history(tail=3)

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.03, stop_time=1.6,
                                  dynamic_substeps=2, rel_fast_threshold=120)

df = dynamics.history.get()
df

# %%
dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['red', 'green', 'gray'])

# %% [markdown]
# ### The (transiently) high value of [U] again led to an increase in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %% [markdown] tags=[]
# # 4. Now, instead, let's DECREASE [U]

# %%
dynamics.set_chem_conc(species_name="U", conc=5., snapshot=True)
dynamics.describe_state()

# %%
dynamics.get_history(tail=3)

# %% [markdown]
# ### Take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.03, stop_time=2.3,
                                  dynamic_substeps=2, rel_fast_threshold=80)

df = dynamics.history.get()
df

# %%
dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['red', 'green', 'gray'])

# %% [markdown]
# ### The (transiently) LOW value of [U] led to an DECREASE in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %%
