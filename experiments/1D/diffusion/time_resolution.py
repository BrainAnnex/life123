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
# **Exploring effect of Time Resolution on accuracy**
#
# In the examples below, the _time advance_ always remains constant,
# but the _number of steps_ used to arrive there vary
#
# NO log file.
#
# LAST REVISED: June 13, 2022

# %%
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home from the folder containing this notebook

# %%
from modules.chemicals.chemicals import Chemicals as chem
from life_1D.bio_sim_1d import BioSim1D as bio


# %%
def set_initial_condition():
    # Set or reset the initial concentrations
    bio.set_uniform_concentration(species_index=0, conc=0.)
    bio.inject_conc_to_bin(bin=2, species_index=0, delta_conc=10.)


# %%
chem_data = chem(diffusion_rates=[0.1])
bio.initialize_system(n_bins=10, chem_data=chem_data)

set_initial_condition()

bio.describe_state(show_diffusion_rates=True)

# %%
t_final = 33.3

bio.diffuse(total_duration=t_final, n_steps=10, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=20, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=30, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=50, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=100, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=1000, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=10000, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=100000, verbose=True)