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
# ### `2 S <-> U` and `S <-> X` (both mostly forward)
#
# ## Based on the experiment `up_regulate_3"
#
# 1st-order kinetics throughout.   
#
# LAST REVISED: Feb. 8, 2023

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

# %%
# Initialize the system
chem_data = chem(names=["U", "X", "S"])

# Reaction 2 S <-> U , with 1st-order kinetics for all species
chem_data.add_reaction(reactants=[(2, "S")], products="U",
                       forward_rate=8., reverse_rate=2.)

# Reaction S <-> X , with 1st-order kinetics for all species
chem_data.add_reaction(reactants="S", products="X",
                       forward_rate=6., reverse_rate=3.)

chem_data.describe_reactions()

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.1, stop_time=0.8)

dynamics.explain_time_advance()

# %%
df = dynamics.get_history()
df

# %%
dynamics.plot_curves(colors=['red', 'green', 'gray'])

# %%
dynamics.diagnostic_data_baselines.get()

# %%
dynamics.get_diagnostic_data(0)

# %%
dynamics.get_diagnostic_data(1)

# %%
