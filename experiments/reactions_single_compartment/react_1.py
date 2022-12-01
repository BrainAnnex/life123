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
# ### A simple A <-> B reaction between 2 species,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# Based on the experiment _"1D/reactions/reaction_1"_ ; this is simply the "single-compartment" version of it.
#
# LAST REVISED: Nov. 30, 2022

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

# %% tags=[]
# Initialize the reaction
chem_data = chem(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

print("Number of reactions: ", chem_data.number_of_reactions())

# %%
chem_data.describe_reactions()

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)

# %%
# Initial concentrations of all the chemicals, in index order
dynamics.set_conc([10., 50.])

# First step of reaction
dynamics.single_compartment_react(time_step=0.1, n_steps=1,
                                  snapshots=None)

# %%
# Numerous more steps
dynamics.single_compartment_react(time_step=0.1, n_steps=10,
                                  snapshots=None)

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order reactions),
#  the systems settles in the following equilibrium:
#
# [A] = 23.99316406
#  
# [B] = 36.00683594
#

# %%
dynamics.is_in_equilibrium(rxn_index=0, conc={"A": dynamics.system[0], "B": dynamics.system[1]})

# %%
