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
# ## `U` ("Up-regulator") up-regulates `X` , by sharing an upstream reagent `S` ("Source") across 2 separate reactions:   
# ### `2 S <-> U` and `S <-> X` (both mostly forward)
#
# 1st-order kinetics throughout.   
#
# Invoking [Le Chatelier's principle](https://www.chemguide.co.uk/physical/equilibria/lechatelier.html), it can be seen that, starting from equilibrium, when [U] goes up, so does [S]; and when [S] goes up, so does [X].   
# Conversely, when [U] goes down, so does [S]; and when [S] goes down, so does [X].   
#
# This experiment is a counterpart of `up_regulate_2`, with "upstream" rather than "downstream" reactions.
#
# LAST REVISED: Feb. 19, 2023

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
chem_data = chem(names=["U", "X", "S"])

# Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants=[(2, "S")], products="U",
                       forward_rate=8., reverse_rate=2.)

# Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants="S", products="X",
                       forward_rate=6., reverse_rate=3.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
dynamics.describe_state()

# %% [markdown] tags=[]
# # 1. Take the initial system to equilibrium

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.01, stop_time=1.5,
                                  dynamic_substeps=4, rel_fast_threshold=50)

df = dynamics.get_history()
#df
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### Note that [S] is initially 0, and that it builds up thru _reverse_ reactions

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=3)

# %% [markdown] tags=[]
# # 2. Now, let's suddenly increase [U]

# %%
dynamics.describe_state()

# %%
dynamics.set_chem_conc(species_name="U", conc=100.)
dynamics.describe_state()

# %%
dynamics.get_history(tail=3)

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.01, stop_time=3.0,
                                  dynamic_substeps=4, rel_fast_threshold=50)

df = dynamics.get_history()
#df
dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### The (transiently) high value of [U] led to an increase in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# # 3. Let's again suddenly increase [U]

# %%
dynamics.describe_state()

# %%
dynamics.set_chem_conc(species_name="U", conc=150.)
dynamics.describe_state()

# %%
dynamics.get_history(tail=3)

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.01, stop_time=4.5,
                                  dynamic_substeps=4, rel_fast_threshold=50)

df = dynamics.history.get()
#df
dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### The (transiently) high value of [U] again led to an increase in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %% [markdown] tags=[]
# # 4. Now, instead, let's DECREASE [U]

# %%
dynamics.describe_state()

# %%
dynamics.set_chem_conc(species_name="U", conc=80.)
dynamics.describe_state()

# %%
dynamics.get_history(tail=3)

# %% [markdown]
# ### Take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.01, stop_time=6.,
                                  dynamic_substeps=4, rel_fast_threshold=50)

df = dynamics.history.get()
#df
dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### The (transiently) LOW value of [U] led to an DECREASE in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %%
