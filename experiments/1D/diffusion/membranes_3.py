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
# # Membranes, with selective permeability, separate previously fully-mixed chemicals
#
# #### Two chemicals: the membranes are permeable to one, and impermeable to the other one 

# %% [markdown]
# ### TAGS :  "membranes 1D", "basic", "diffusion 1D"

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
# with two chemicals `A` and `B`

# %%
chem_data = ChemData(diffusion_rates=[10., 5.], plot_colors=["turquoise", "green"])   # Name "A", "B" automatically assigned

bio = BioSim1D(n_bins=7, chem_data=chem_data)

# %%
bio.set_bin_conc(bin_address=2, conc=100., chem_label="A") 
bio.set_bin_conc(bin_address=2, conc=100., chem_label="B")

bio.set_membranes(membranes=[ (2,3) ])    # By default impermeable
bio.change_permeability("B", 1.)          # Permeable to `B` (and only to `B`)

bio.describe_state()

# %%
bio.system_heatmaps()

# %%



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
# ### The concentrations of `A` and `B` start to differ
# Membranes are impermeable to `A` but permeable to `B`
#

# %%

# %% [markdown]
# ### Let's advance the diffusion

# %%
bio.diffuse(time_step=0.02, n_steps=4)

# %%
bio.system_heatmaps()

# %%
bio.visualize_system()   # Line curve view

# %%

# %%
bio.diffuse(time_step=0.02, n_steps=15)

# %%
bio.system_heatmaps()

# %%
bio.visualize_system()   # Line curve view

# %%

# %%
bio.diffuse(time_step=0.02, n_steps=700)

# %%
bio.system_heatmaps()

# %%
bio.visualize_system()   # Line curve view

# %% [markdown]
# # The membranes, with their selective permeability, have managed to separate the previously fully-mixed `A` and `B`|

# %%

# %% [markdown]
# #### Verify that the total amount of each chemicals hasn't changed; in the case of `B`, it has simply shifted across bins

# %%
bio.check_mass_conservation(expected=100, chem_label="A")

# %%
bio.check_mass_conservation(expected=100, chem_label="B")

# %%
