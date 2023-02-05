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
# ## Diffusion of a bell-shaped initial concentration along an initial gradient
# A single chemical, whose initial concentration is a mix of a bell shape and a gradient.   
# Contrary to perhaps an intuition of a "pile sliding down a sand dune as a unit", the concentration peak
# remains in place, and simply spreads out from there
#
# LAST REVISED: Nov. 28, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from src.life_1D.bio_sim_1d import BioSim1D

import plotly.express as px
from src.modules.reactions.reaction_data import ReactionData as chem

# %%
# Initialize the system
chem_data = chem(names=["A"], diffusion_rates=[10.])

bio = BioSim1D(n_bins=200, chem_data=chem_data)

# %%
# Set up the initial bell-shape concentration, with the peak close to one end of the system
bio.inject_bell_curve(species_name="A", mean=0.25, sd=0.1, amplitude=20., bias=0)

# %%
# Visualize the system state so far
bio.visualize_system(caption="Preparations in progress")

# %%
# Complete the initial-system preparation by adding a gradient slanting to the right
bio.inject_gradient("A", conc_left = 50., conc_right = 0.)

# %%
# Visualize the complete initial state
bio.visualize_system()

# %%
# Do several round of diffusion, over relatively small times
for _ in range(5):
    bio.diffuse(total_duration=20, n_steps=1000)
    bio.visualize_system()

# %%
# Do more rounds of diffusion, over large times spans
for _ in range(9):
    bio.diffuse(total_duration=300, n_steps=10000)
    bio.visualize_system()

# %% [markdown]
# The system is **approaching equilibrium** at a value a little less than 45.   
#
# That's consisten with the contribution from:  
# 1) a slightly truncated bell curve (whose total un-truncated area underneath is 1, multiplied by the "amplitude" factor of 20)
# 2) the 50-0 gradient (with average value 25)

# %%
