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
# LAST REVISED: Aug. 30, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home from the folder containing this notebook

# %%
from modules.chemicals.chemicals import Chemicals as chem
from life_1D.bio_sim_1d import BioSim1D as bio

# %%
chem_data = chem(names=["A", "B", "C"])
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
