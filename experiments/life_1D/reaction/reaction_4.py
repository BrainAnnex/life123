# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# **One-bin A + B <-> C, with 1st-order kinetics for each species,
# taken to equilibrium**
#
# Diffusion not applicable (just 1 bin)
#
# LAST REVISED: July 2, 2022

# %%
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "C"])     # NOTE: Diffusion not applicable (just 1 bin)

rxn = Reactions(chem_data)

# Reaction A + B <-> C , with 1st-order kinetics for each species
rxn.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                 forward_rate=5., reverse_rate=2.)

bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)
bio.set_uniform_concentration(species_index=2, conc=20.)

bio.describe_state()

# %%
rxn.describe_reactions()

# %%
rxn._internal_reactions_data()    # Low-level view of the reactions data

# %%
# First step
bio.react(time_step=0.002, n_steps=1)
bio.describe_state()

# %% [markdown]
# _Early in the reaction :_
# [A] = 5.08 ,  [B] = 45.08 ,  [C] = [24.92]

# %%
# Numerous more steps
bio.react(time_step=0.002, n_steps=40)

bio.describe_state()

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  
# [A] = 0.29487741 , [B] = 40.29487741 , [C] = [29.70512259]

# %% [markdown]
# # Note: "A" (now almost completely depleted) is largely the limiting reagent

# %%
A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
print(f"Ratio of equilibrium concentrations (C_eq / (A_eq * B_eq)) : {C_eq / (A_eq * B_eq)}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")

# %%
