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
# ## Exploring reaching equilibrium, first on a shorter timescale and then a longer one 
# ### (but both with identical time steps.)
#
# The system starts out with a "concentration pulse" in bin 2 (the 3rd bin from the left) - i.e. that bin is initially the only one with a non-zero concentration of the only chemical species.
#
# Notice the diffusing pulse "bouncing" off the left wall after total time 30.
#
# Then the system is left undisturbed, and followed to equilibrium.
#
# *NO log file.*
#
# LAST REVISED: June 4, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from src.modules.reactions.reaction_data import ChemData as chem
from src.life_1D.bio_sim_1d import BioSim1D

# %%
chem_data = chem(diffusion_rates=[0.1])
bio = BioSim1D(n_bins=10, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=0.)
bio.inject_conc_to_bin(species_index=0, bin_address=2, delta_conc=10.)

bio.describe_state()

# %%
print("\n\nSTARTING on SHORTER time scales.  Dtime=10, with time steps of 0.1 ...")

# %%
for i in range(10):
    delta_time = 10.
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)
    print(f"\nAfter Delta time {delta_time}.  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)

# %% [markdown]
# *Notice the diffusing pulse "bounces" off the left wall after total time 30:*  
# the concentration at cell 0 increases from t=0 to 30, and then it's coming down by t=40

# %%
print("\n\nREPEATING to LONGER time scales.  Dtime=100, again with time steps of 0.1 ...")

# Reset the concentrations
bio.set_uniform_concentration(species_index=0, conc=0.)
bio.inject_conc_to_bin(species_index=0, bin_address=2, delta_conc=10.)

#total_time = 0.
for i in range(20):
    delta_time = 100.
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)
    print(f"\nAfter Delta time {delta_time}.  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)

# %% [markdown]
# ## The system has now reached equilibrium

# %%
