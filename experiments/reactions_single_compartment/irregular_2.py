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
# # Continuation of "irregular_1"
#
# ## A cycle of reactions A <-> B <-> C <-> A, 
# The "closing" of the cycle (the "return" parth from C to A) is couple with an "energy donor" reaction:
# C + E_High <-> A + E_Low
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

# %%

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

# %%
initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}
initial_conc

# %% [markdown]
# # START OVER WITH COARSE FIXED-RESOLUTION

# %%
dynamics2 = ReactionDynamics(reaction_data=chem_data)
dynamics2.set_conc(conc= initial_conc,
                  snapshot=True)      # Note the abundant energy source
dynamics2.describe_state()

# %%
dynamics2.single_compartment_react(time_step=0.0005,stop_time=0.005)
dynamics2.single_compartment_react(time_step=0.001,stop_time=0.2)

# %%
dynamics2.single_compartment_react(time_step=0.002, stop_time=3.)

# %%
dynamics2.single_compartment_react(time_step=0.005, stop_time=8.)


# %%
def graphics():
    df = dynamics2.history.get()
    fig = px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B", "C"], 
                  title="Changes in concentrations",
                  color_discrete_sequence = ['blue', 'green', 'brown'],
                  labels={"value":"concentration", "variable":"Chemical"})
    fig.show()   


# %%
graphics()

# %%
