# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Membranes in 1D : Diffusion and Passive Transport across Membranes
#
# #### Simple scenarios

# %% [markdown]
# ### TAGS :  "membranes 1D", "basic", "quick-start"

# %%
LAST_REVISED = "May 18, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from life123 import BioSim1D, ChemData, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Prepare the initial system
# with a single chemical `A`

# %%
chem_data = ChemData(diffusion_rates=10.)   # Name "A" automatically assigned to single chemical

bio = BioSim1D(n_bins=4, chem_data=chem_data)

# %%
initial_conc = [20, 50, 150, 250]

# %%
bio.set_species_conc(conc_list=initial_conc, chem_index=0) 

bio.describe_state()

# %%
bio.system_snapshot()

# %%
bio.visualize_system()   # Line curve view

# %%
bio.system_heatmaps()

# %%



# %% [markdown]
# # Add Membranes

# %% [markdown]
# #### The bins would normally mix up by diffusion, but we'll prevent that by inserting impermeable membranes

# %%
bio.set_membranes(membranes=[ (1,2) ])    # By default impermeable

# %%
bio.membranes

# %%
bio.describe_state()

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ## Now, let's start the diffusion

# %%
bio.diffuse(time_step=0.02, n_steps=4)

# %%
bio.system_snapshot()

# %% [markdown]
# ### Nothing has changed in the 2 leftmost bins, as a result of the impermeable membranes  
# By contrast, diffusion is progressing between the 2 rightmost bins
#

# %%
bio.system_heatmaps()

# %%
bio.diffuse(time_step=0.02, n_steps=16)

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ## Now, let's assign a small permeability to our single chemical `A` to the membranes (by passive transport)

# %%
bio.change_permeability("A", 1.)

# %%
bio.diffuse(time_step=0.02, n_steps=1)

# %%
bio.system_heatmaps()

# %% [markdown]
# ### Passive transport across membranes is now taking place both out of bin 1 (to its left neighbor) and into it (from its right neighbors)  
# Diffusion is continuing normally between the 2 rightmost bins

# %%
bio.diffuse(time_step=0.02, n_steps=4)

# %%
bio.system_heatmaps()

# %%
bio.diffuse(time_step=0.02, n_steps=45)

# %%
bio.system_heatmaps()

# %%
bio.diffuse(time_step=0.02, n_steps=500)

# %%
bio.system_heatmaps()

# %% [markdown]
# #### Verify that the total amount of our chemical present hasn't changed; it has simply shifted across bins

# %%
bio.check_mass_conservation(expected=sum(initial_conc), chem_index=0)

# %%
