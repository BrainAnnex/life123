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
# ## A cycle of reactions A <-> B <-> C <-> A, 
# ### The "closing" of the cycle (the "return" parth from C to A) is couple with an "energy donor" reaction:
# ### C + E_High <-> A + E_Low
# where E_High and E_Low are, respectively, the high- and low- energy molecules that drive the cycle (for example, think of ATP/ADP)
#
# All 1st-order kinetics.   
#
# LAST REVISED: Jan. 21, 2023

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
chem_data = chem(names=["A", "B", "C", "E_high", "E_low"])

# Reaction A <-> B, mostly in forward direction (favored energetically)
# Note: all reactions in this experiment have 1st-order kinetics for all species
chem_data.add_reaction(reactants="A", products="B",
                       forward_rate=9., reverse_rate=3.)

# Reaction B <-> C, also favored energetically
chem_data.add_reaction(reactants="B", products="C",
                       forward_rate=8., reverse_rate=4.)

# Reaction C + E_High <-> A + E_Low, also favored energetically, but kinetically slow
# Note that, thanks to the energy donation from E, we can go "upstream" from C, to the higher-energy level of "A"
chem_data.add_reaction(reactants=["C" , "E_high"], products=["A", "E_low"],
                       forward_rate=1., reverse_rate=0.2)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}
initial_conc

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc=initial_conc,
                  snapshot=True)      # Note the abundant energy source
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Start the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.001, n_steps=30, 
                                  dynamic_steps=4, fast_threshold=10)

# %%
dynamics.history.get()

# %%
dynamics.explain_time_advance()

# %%
dynamics.single_compartment_react(time_step=0.01, stop_time=8.,
                                  dynamic_steps=8, fast_threshold=10)

df = dynamics.history.get()
df

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['blue', 'green', 'brown'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
fig = px.line(data_frame=df, x="SYSTEM TIME", y=["E_high", "E_low"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['red', 'gray'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
dynamics.is_in_equilibrium(tolerance=4)

# %%
df.head(22)

# %%
df[113:131]

# %%
df[df['SYSTEM TIME'] > 2.52].head(50)

# %%
dynamics.diagnostic_data_baselines.get()[100:130]

# %%
