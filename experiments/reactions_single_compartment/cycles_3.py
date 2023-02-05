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
# # IN-PROGRESS
# ## A cycle of reactions `A <-> B <-> C <-> A`
# #### the "closing" of the above cycle (the "return" parth from `C` to `A`) is coupled with an "energy donor" reaction:
# #### `C + E_High <-> A + E_Low`
# #### where `E_High` and `E_Low` are, respectively, the high- and low- energy molecules that drive the cycle (for example, think of ATP/ADP) 
# ### AND `E_High` is periodically replenished
#
# All 1st-order kinetics.    
#
# LAST REVISED: Feb. 1, 2023

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.modules.numerical.numerical import Numerical as num

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
                       forward_rate=1., reverse_rate=0.5)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}  # Note the abundant energy source "E_high"
initial_conc

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc=initial_conc, snapshot=True)
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %%
dynamics.single_compartment_react(time_step=0.0016, reaction_duration=0.025,
                                  dynamic_substeps=2, rel_fast_threshold=150.)
#dynamics.get_history()

# %%
dynamics.describe_state()

# %%
dynamics.explain_time_advance()

# %%
dynamics.set_chem_conc(conc=1000., species_name="E_high", snapshot=True)
dynamics.set_chem_conc(conc=0., species_name="E_low", snapshot=True)

# %%
for i in range(19):
    dynamics.single_compartment_react(time_step=0.001, reaction_duration=0.025,
                                      dynamic_substeps=2, rel_fast_threshold=150., silent=True)
    dynamics.set_chem_conc(conc=1000., species_name="E_high", snapshot=True)
    dynamics.set_chem_conc(conc=0., species_name="E_low", snapshot=True)

# %%
dynamics.system_time

# %%
for i in range(30):
    dynamics.single_compartment_react(time_step=0.001, reaction_duration=0.025,
                                      dynamic_substeps=2, rel_fast_threshold=150., silent=True)
    dynamics.set_chem_conc(conc=1000., species_name="E_high", snapshot=True)
    dynamics.set_chem_conc(conc=0., species_name="E_low", snapshot=True)

# %%
dynamics.system_time

# %%
for i in range(250):
    dynamics.single_compartment_react(time_step=0.001, reaction_duration=0.025,
                                      dynamic_substeps=2, rel_fast_threshold=120., silent=True)
    dynamics.set_chem_conc(conc=1000., species_name="E_high", snapshot=True)
    dynamics.set_chem_conc(conc=0., species_name="E_low", snapshot=True)
    if i % 10 == 0:
        print(f"Processed {i} rounds so far...")

# %%
dynamics.plot_curves(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_curves(chemicals=["A", "B", "C"])

# %%
dynamics.is_in_equilibrium()

# %%
