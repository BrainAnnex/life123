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
# **THE DANGER OF EXCESSIVELY LARGE SINGLE TIME STEPS IN DIFFUSION**
#
# (We're assuming bins of size 1 ; TODO: adjust for any size)
#
# ->  When the time step is (0.5 / diffusion rate),
#     a 2-bin system equilibrates in a single step,
#     and some 3-bin systems can over-shoot equilibrium!
#
# ->  When the time step is (0.33333 / diffusion rate),
#     some 3-bin systems equilibrate in a single step
#
# So, (0.33333 / diffusion rate) is a - rather lax - upper bound for
# sensible single time steps!
#
# That value of 0.33333 is saved in the Class variable "time_step_threshold"
#
# *NO log file.*
#
# LAST REVISED: Aug. 28, 2022

# %%
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home from the folder containing this notebook

# %%
from modules.chemicals.chemicals import Chemicals as chem
from life_1D.bio_sim_1d import BioSim1D as bio

# %%
###########################
# Simulate a 2-bin system
###########################
chem_data = chem(diffusion_rates=[10.])
bio.initialize_system(n_bins=2, chem_data=chem_data)
bio.inject_conc_to_bin(bin_address=0, delta_conc=100., species_index=0)
bio.describe_state()

# %%
bio.time_step_threshold = 0.51  # To bypass limits that are typically in place

# When the time step is (0.5 / diffusion rate),
# a 2-bin system equilibrates in a single step!
bio.diffuse(time_step=0.05, n_steps=1)
print(bio.system) 

# %% [markdown]
# _Note: [[50. 50.]]  : **the two bins have equilibrated!**_

# %%
######################################
# Simulate a 3-bin system
# with an excessive single time step
######################################
bio.initialize_system(n_bins=3, chem_data=chem_data)
bio.inject_conc_to_bin(bin_address=1, delta_conc=100., species_index=0)
bio.describe_state()
#3 bins and 1 species:   [[  0. 100.   0.]]

# %%
bio.time_step_threshold = 0.51      # To bypass a limit that is typically in place

# When the time step is (0.5 / diffusion rate),
# a 3-bin system can overshoot equilibrium!
bio.diffuse(time_step=0.05, n_steps=1)
print(bio.system)

# %% [markdown]
# _Note: [[50.  0. 50.]] : the diffusion has **over-shot equilibrium!!!**_

# %%
#############################################
# Re-Simulate the 3-bin system,
# with a somewhat smaller single time step
#############################################
bio.initialize_system(n_bins=3, chem_data=chem_data)
bio.inject_conc_to_bin(bin_address=1, delta_conc=100., species_index=0)
bio.describe_state()

# %%
bio.time_step_threshold = 0.34      # To bypass a limit that is typically in place

# When the time step is (0.33333 / diffusion rate),
# a 3-bin system, configured as above, equilibrates in a single step!
bio.diffuse(time_step=0.033333, n_steps=1)
print(bio.system)  

# %% [markdown]
# _Note: the three bins have **equilibrated**!_
