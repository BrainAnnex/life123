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
# #### Simple scenarios with 1 chemical

# %% [markdown]
# ### TAGS :  "membranes 1D", "basic", "quick-start", "diffusion 1D"

# %%
LAST_REVISED = "May 19, 2025"
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
bio.diffuse(time_step=0.02, n_steps=1)

# %%
bio.system_snapshot()

# %%
bio.system_heatmaps()

# %% [markdown]
# ### Nothing has changed in the 2 leftmost bins, as a result of the impermeable membranes  
# By contrast, diffusion is progressing between the 2 rightmost bins
#

# %% [markdown]
# *Technical side note:* since we're using, by default, simple 
# 3-1 stencils for the diffusion step, the concentration increment in bin 3, and its corresponding decrement in bin 2, is:

# %%
0.02 * 100 * 10 / (1*1)    # (Time step) * (Delta concentration) * (Diffusion rate const) / (bin length squared)

# %% [markdown]
# The 150 initial value in bin 2 increased by 20 to 170, and the 250 value in bin 3 correspondingly decreased to 230

# %%

# %%

# %% [markdown]
# ### Let's further advance the diffusion, by many more steps

# %%
bio.diffuse(time_step=0.02, n_steps=14)

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

# %% [markdown]
# *Technical side note -* Passive membrane transport out of bin 1 (conc. 50) into bin 0 (conc. 20) in this single step was:

# %%
30 * 1 * 0.02   # (Delta conc) * Permeability * (Time step) , since all bin sizes are 1

# %% [markdown]
# The concentration of bin 0 has increased from 20 to 20.6 during this single step.  
# The concentration of bin 1 is also increasing, because its loss to bin 0 is more than compensated by its gain from bin 2 (bigger conc. difference)

# %%

# %% [markdown]
# ### Let's continue to advance the diffusion

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
