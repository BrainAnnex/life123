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
# This experiment is a counterpart of experiment `up_regulate_2`, with "upstream" rather than "downstream" reactions.
#
# Note: numerical errors in the same reactions (with the same initial conditions) is explored in the experiment "large_time_steps_2"
#
# LAST REVISED: May 26, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ChemData as chem
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

# %%

# %% [markdown] tags=[]
# # 1. Take the initial system to equilibrium

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.5, downshift=0.5, abort=0.5)
dynamics.set_error_step_factor(0.5)

dynamics.single_compartment_react(initial_step=0.01, target_end_time=1.5,
                                  variable_steps=True, explain_variable_steps=False)

#df = dynamics.get_history()
#dynamics.explain_time_advance()

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
dynamics.is_in_equilibrium()

# %%

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
dynamics.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4)  # Needs to tighten to time advance, to prevent mild instability


dynamics.single_compartment_react(initial_step=0.01, target_end_time=3.0,
                                  variable_steps=True, explain_variable_steps=False)

#df = dynamics.get_history()
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### The (transiently) high value of [U] led to an increase in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%

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
dynamics.single_compartment_react(initial_step=0.01, target_end_time=4.5,
                                  variable_steps=True, explain_variable_steps=False)

#dynamics.get_history()

#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### The (transiently) high value of [U] again led to an increase in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %%

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
dynamics.single_compartment_react(initial_step=0.01, target_end_time=6.,
                                  variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.history.get_dataframe()
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### The (transiently) LOW value of [U] led to an DECREASE in [X]

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %% [markdown]
# **IDEAS TO EXPLORE**:   
#
# * Effect of the stoichiometry and the Delta_G on the "amplification" of the signal (from [U] to [X]) 
#
# * Effect of a continuously-varying (maybe oscillating [U]), and its being affected by the reactions' kinetics
#
# * Combining this experiment and `up_regulate_2` in a *"bifan motif"*

# %%
