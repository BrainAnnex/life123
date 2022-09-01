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
# ### Exploring the data structures of MEMBRANES, and reactions in them 
# #### - with NO DIFFUSION
#
# LAST REVISED: Aug. 31, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home from the folder containing this notebook

# %%
from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions

from life_1D.bio_sim_1d import BioSim1D as bio

# %%
chem_data = chem(names=["A", "B", "C"])     # NOTE: Diffusion not done
bio.initialize_system(n_bins=5, chem_data=chem_data)

bio.set_membranes(membrane_pos=[1])   # A single membrane, passing thru bin 1

bio.set_all_uniform_concentrations(conc_list=[4., 8., 12.])

bio.describe_state()

# %%
bio.set_bin_conc(bin_address=1, species_name="A", conc=10.)
bio.set_bin_conc(bin_address=1, species_name="A", conc=55., across_membrane=True)
bio.set_bin_conc(bin_address=1, species_name="B", conc=20.)
bio.set_bin_conc(bin_address=1, species_name="C", conc=30., both_sides=True)

bio.describe_state()

# %%
# Make the last bin match all the concentrations of the "post-membrane" section of bin 1
bio.set_bin_conc(bin_address=4, species_name="A", conc=55.)
bio.set_bin_conc(bin_address=4, species_name="C", conc=30.)

bio.describe_state()

# %%
# Reaction A + B <-> C , with 1st-order kinetics in both directions, mostly forward
bio.all_reactions.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=8., reverse_rate=2.)

bio.all_reactions.describe_reactions()

# %%
bio.react(time_step=0.002, n_steps=1)
bio.describe_state()

# %% [markdown]
# ### Note how (in the absence of diffusion, which we're neglecting) the concentrations on the "post-membrane side" of bin 1 continue to match those of bin 5

# %% [markdown]
# ## Now continue to reaction equilibrium

# %%
bio.react(time_step=0.002, n_steps=100)
bio.describe_state()

# %% [markdown]
# ### The system has now reached equilibrium
# ### in individual bins, which remain separate because we're NOT doing diffusion in this experiment

# %% [markdown]
# Verify the equilibrium in each of the active bins

# %%
bio.all_reactions.is_in_equilibrium(rxn_index=0, conc={"A": 0.7932534523016195, "B": 4.793253452301622, "C": 15.206746547698396})
# A was largely the limiting reagent

# %%
bio.all_reactions.is_in_equilibrium(rxn_index=0, conc={"A": 0.897094737056602, "B": 10.897094737056594, "C": 39.10290526294337})
# A was largely the limiting reagent

# %%
bio.all_reactions.is_in_equilibrium(rxn_index=0, conc={"A": 47.20020986266437, "B": 0.20020986266437912, "C": 37.79979013733563})
# This time, with ample [A], the limiting reagent was B

# %%
