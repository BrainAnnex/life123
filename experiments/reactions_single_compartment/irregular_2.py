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
initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}  # Note the abundant energy source
initial_conc

# %% [markdown]
# # _COARSE_ FIXED-RESOLUTION

# %%
dynamics2 = ReactionDynamics(reaction_data=chem_data)
dynamics2.set_conc(conc= initial_conc,
                  snapshot=True)      
dynamics2.describe_state()

# %%
# time_step values are picked as large as possible, but short of causing errors
dynamics2.single_compartment_react(time_step=0.0005,stop_time=0.005)

# %%
dynamics2.single_compartment_react(time_step=0.001,stop_time=0.2)

# %%
dynamics2.single_compartment_react(time_step=0.002, stop_time=3.)

# %%
dynamics2.single_compartment_react(time_step=0.005, stop_time=8.)

# %%
dynamics2.plot_curves()

# %%
df2 = dynamics2.get_history()
len(df2)

# %%
df2["A"] > df2["B"]


# %%
def foo(df, t_min, t_max, var1, var2):
    df


# %%
(df2["SYSTEM TIME"] > 2.) & (df2["SYSTEM TIME"] < 2.4)

# %%
df2[df2.eval("`SYSTEM TIME` >= 2.31 and `SYSTEM TIME` <= 2.33")]

# %%
df2[df2.eval("`SYSTEM TIME` >= 2.31 and `SYSTEM TIME` <= 2.33")]

# %%
df2[df2["SYSTEM TIME"].between(2.31, 2.33)]

# %%
df2[(df2["SYSTEM TIME"] >= 2.31) & (df2["SYSTEM TIME"] <= 2.33)]

# %% [markdown]
# # START OVER WITH _FINE_ FIXED-RESOLUTION

# %%
dynamics3 = ReactionDynamics(reaction_data=chem_data)
dynamics3.set_conc(conc= initial_conc,
                  snapshot=True)      
dynamics3.describe_state()

# %%
# Using the smallest time_step value that was ever used in the subset of the VARIABLE-resolution run
dynamics3.single_compartment_react(time_step=0.00025, stop_time=8.)

# %%
dynamics3.plot_curves()

# %%
len(dynamics3.history.get())

# %%
